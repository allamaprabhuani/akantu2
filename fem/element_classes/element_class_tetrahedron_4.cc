/**
 * @file   element_class_tetrahedron_4.cc
 * @author Guillaume ANCIAUX <anciaux@epfl.ch>
 * @date   Mon Aug 16 18:09:53 2010
 *
 * @brief  Specialization of the element_class class for the type _tetrahedron_4
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
     x (0,0,1,0)
     |`
     |  `  °             \xi
     |    `    °       -
     |      `       x (0,0,0,1)
     |      q.`   -  '
     |         -`     '
     |     -       `   '
     | -             `  '
     x------------------x-----> \zeta
 (1,0,0,0)           (0,1,0,0)
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta - \zeta \\
 * N2 &=& \xi \\
 * N3 &=& \eta \\
 * N4 &=& \zeta
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f[
 * \xi_{q0} = 1/4 \qquad  \eta_{q0} = 1/4  \qquad  \zeta_{q0} = 1/4
 * @f]
 */

/* -------------------------------------------------------------------------- */

  // /// shape functions
  // shape[0] = 1./4.; /// N1(q_0)
  // shape[1] = 1./4.; /// N2(q_0)
  // shape[2] = 1./4.; /// N3(q_0)
  // shape[3] = 1./4.; /// N4(q_0)


/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedron_4>::nb_nodes_per_element;
template<> UInt ElementClass<_tetrahedron_4>::nb_quadrature_points;
template<> UInt ElementClass<_tetrahedron_4>::spatial_dimension;


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_tetrahedron_4>::computeShapes(const Real * natural_coords, 
								   Real * shapes){
  Real c0 = natural_coords[1]; /// @f$ c0 = \eta @f$
  Real c1 = natural_coords[2]; /// @f$ c1 = \zeta @f$
  Real c2 = 1 - natural_coords[0] -  natural_coords[1] -  natural_coords[2];/// @f$ c2 = 1 - \xi - \eta - \zeta @f$
  Real c3 = natural_coords[0]; /// @f$ c2 = \xi @f$
  
  shapes[0] = c0;
  shapes[1] = c1;
  shapes[2] = c2;
  shapes[3] = c3;
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_tetrahedron_4>::computeDNDS(__attribute__ ((unused)) const Real * natural_coords,
								 Real * dnds) {

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccc}
   *            \frac{\partial N1}{\partial \xi}   & \frac{\partial N2}{\partial \xi}
   *          & \frac{\partial N3}{\partial \xi}   & \frac{\partial N4}{\partial \xi} \\
   *            \frac{\partial N1}{\partial \eta}  & \frac{\partial N2}{\partial \eta}
   *          & \frac{\partial N3}{\partial \eta}  & \frac{\partial N4}{\partial \eta} \\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial \zeta}
   *          & \frac{\partial N3}{\partial \zeta} & \frac{\partial N4}{\partial \zeta}
   *          \end{array}
   *        \right)
   * @f]
   */

  dnds[0] = -1.; dnds[1] = 1.; dnds[2]  = 0.; dnds[3]  = 0.;
  dnds[4] = -1.; dnds[5] = 0.; dnds[6]  = 1.; dnds[7]  = 0.;
  dnds[8] = -1.; dnds[9] = 0.; dnds[10] = 0.; dnds[11] = 1.;


}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_tetrahedron_4>::computeJacobian(const Real * dxds,
								     const UInt dimension, 
								     Real & jac) {

  if (dimension == spatial_dimension){
    Real weight = 1./6.;
    Real det_dxds = Math::det3(dxds);
    jac = det_dxds * weight;    
  }  
  else {
    AKANTU_DEBUG_ERROR("to be implemented");
  }
}
 
/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_tetrahedron_4>::getInradius(const Real * coord) {
  return Math::tetrahedron_inradius(coord, coord+3, coord+6, coord+9);
}
