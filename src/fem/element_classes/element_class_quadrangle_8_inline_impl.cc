/**
 * @file   element_class_quadrangle_8_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue May 17 18:35:04 2011
 *
 * @brief  Specialization of the ElementClass for the _quadrangle_8
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
	       |
   (-1,1)    (0,1)   (1,1)
       x-------x-------x
       |       |       |
       |       |       |
       |       |       |
 (-1,0)|       |       |(1,0)
   ----x---------------X----->  \xi
       |       |       |
       |       |       |
       |       |       |
       |       |       |
       x-------x-------x
   (-1,-1)   (0,-1)  (1,-1)
	       |
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{lll}
 * N1 = (1 - \xi) (1 - \eta)(- 1 - \xi - \eta) / 4
 *       & \frac{\partial N1}{\partial \xi}  = (1 - \eta)(2 \xi + \eta) / 4
 *       & \frac{\partial N1}{\partial \eta} = (1 - \xi)(\xi + 2 \eta) / 4 \\
 * N2 = (1 + \xi) (1 - \eta)(- 1 + \xi - \eta) / 4 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta)(2 \xi - \eta) / 4
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi)(\xi - 2 \eta) / 4 \\
 * N3 = (1 + \xi) (1 + \eta)(- 1 + \xi + \eta) / 4 \\
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta)(2 \xi + \eta) / 4
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi)(\xi + 2 \eta) / 4 \\
 * N4 = (1 - \xi) (1 + \eta)(- 1 - \xi + \eta) / 4
 *       & \frac{\partial N4}{\partial \xi}  = (1 + \eta)(2 \xi - \eta) / 4
 *       & \frac{\partial N4}{\partial \eta} = - (1 - \xi)(\xi - 2 \eta) / 4 \\
 * N5 = (1 - \xi^2) (1 - \eta) / 2
 *       & \frac{\partial N1}{\partial \xi}  = - \xi (1 - \eta)
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi^2) / 2  \\
 * N6 = (1 + \xi) (1 - \eta^2) / 2 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta^2) / 2
 *       & \frac{\partial N2}{\partial \eta} = - \eta (1 + \xi) \\
 * N7 = (1 - \xi^2) (1 + \eta) / 2 \\
 *       & \frac{\partial N3}{\partial \xi}  = - \xi (1 + \eta)
 *       & \frac{\partial N3}{\partial \eta} = (1 - \xi^2) / 2 \\
 * N8 = (1 - \xi) (1 - \eta^2) / 2
 *       & \frac{\partial N4}{\partial \xi}  = - (1 - \eta^2) / 2
 *       & \frac{\partial N4}{\partial \eta} = - \eta (1 - \xi) \\
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 0 \qquad  \eta_{q0} = 0
 * @f}
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_quadrangle_8>::nb_nodes_per_element;
template<> UInt ElementClass<_quadrangle_8>::nb_quadrature_points;
template<> UInt ElementClass<_quadrangle_8>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_quadrangle_8>::computeShapes(const Real * natural_coords,
								   Real * shapes) {
  /// Natural coordinates
  const Real xi  = natural_coords[0];
  const Real eta = natural_coords[1];

  shapes[0] = .25 * (1 - xi) * (1 - eta) * (- 1 - xi - eta);
  shapes[1] = .25 * (1 + xi) * (1 - eta) * (- 1 + xi - eta);
  shapes[2] = .25 * (1 + xi) * (1 + eta) * (- 1 + xi + eta);
  shapes[3] = .25 * (1 - xi) * (1 + eta) * (- 1 - xi + eta);
  shapes[4] =  .5 * (1 - xi * xi) * (1 - eta      );
  shapes[5] =  .5 * (1 + xi     ) * (1 - eta * eta);
  shapes[6] =  .5 * (1 - xi * xi) * (1 + eta      );
  shapes[7] =  .5 * (1 - xi     ) * (1 - eta * eta);
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_quadrangle_8>::computeDNDS(const Real * natural_coords,
								 Real * dnds) {
  const Real xi  = natural_coords[0];
  const Real eta = natural_coords[1];

  /// dN/dxi
  dnds[0]  = .25 * (1 - eta) * (2 * xi + eta);
  dnds[1]  = .25 * (1 - eta) * (2 * xi - eta);
  dnds[2]  = .25 * (1 + eta) * (2 * xi + eta);
  dnds[3]  = .25 * (1 + eta) * (2 * xi - eta);
  dnds[4]  = - xi * (1 - eta);
  dnds[5]  =   .5 * (1 - eta * eta);
  dnds[6]  = - xi * (1 + eta);
  dnds[7]  = - .5 * (1 - eta * eta);

  /// dN/deta
  dnds[8]  = .25 * (1 - xi) * (2 * eta + xi);
  dnds[9]  = .25 * (1 + xi) * (2 * eta - xi);
  dnds[10] = .25 * (1 + xi) * (2 * eta + xi);
  dnds[11] = .25 * (1 - xi) * (2 * eta - xi);
  dnds[12] = -  .5 * (1 - xi * xi);
  dnds[13] = - eta * (1 + xi);
  dnds[14] =    .5 * (1 - xi * xi);
  dnds[15] = - eta * (1 - xi);
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_quadrangle_8>::computeJacobian(const Real * dxds,
								     const UInt dimension,
								     Real & jac){
  if (dimension == spatial_dimension){
    Real det_dxds = Math::det2(dxds);
    jac = det_dxds;
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}



/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_quadrangle_8>::getInradius(const Real * coord) {
  Real a, b, h;

  a = Math::distance_2d(coord + 0*2, coord + 4*2);
  b = Math::distance_2d(coord + 4*2, coord + 1*2);
  h = std::min(a, b);

  a = Math::distance_2d(coord + 1*2, coord + 5*2);
  b = Math::distance_2d(coord + 5*2, coord + 2*2);
  h = std::min(h, std::min(a, b));

  a = Math::distance_2d(coord + 2*2, coord + 6*2);
  b = Math::distance_2d(coord + 6*2, coord + 3*2);
  h = std::min(h, std::min(a, b));

  a = Math::distance_2d(coord + 3*2, coord + 7*2);
  b = Math::distance_2d(coord + 7*2, coord + 0*2);
  h = std::min(h, std::min(a, b));


  return h;
}
