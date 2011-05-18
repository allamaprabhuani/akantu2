/**
 * @file   element_class_quadrangle_8.cc
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
 * @f{equerry*}{
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
  const Real * c = natural_coords;

  shapes[0] = - .25 * (1 - c[0]) * (1 - c[1]) * (1 + c[0] + c[1]);
  shapes[1] = - .25 * (1 + c[0]) * (1 - c[1]) * (1 - c[0] + c[1]);
  shapes[2] = - .25 * (1 + c[0]) * (1 + c[1]) * (1 - c[0] - c[1]);
  shapes[3] = - .25 * (1 - c[0]) * (1 + c[1]) * (1 + c[0] - c[1]);
  shapes[4] =    .5 * (1 - c[0] * c[0]) * (1 - c[1]       );
  shapes[5] =    .5 * (1 + c[0]       ) * (1 - c[1] * c[1]);
  shapes[6] =    .5 * (1 - c[0] * c[0]) * (1 + c[1]       );
  shapes[7] =    .5 * (1 - c[0]       ) * (1 - c[1] * c[1]);
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_quadrangle_8>::computeDNDS(const Real * natural_coords,
								 Real * dnds) {
  const Real * c = natural_coords;

  dnds[0]  = .25 * (1 - c[1]) * (2 * c[0] + c[1]);
  dnds[1]  = .25 * (1 - c[1]) * (2 * c[0] - c[1]);
  dnds[2]  = .25 * (1 + c[1]) * (2 * c[0] + c[1]);
  dnds[3]  = .25 * (1 + c[1]) * (2 * c[0] - c[1]);
  dnds[4]  = - c[0] * (1 - c[1]);
  dnds[5]  =     .5 * (1 - c[1] * c[1]);
  dnds[6]  = - c[0] * (1 + c[1]);
  dnds[7]  = -   .5 * (1 - c[1] * c[1]);

  dnds[8]  = .25 * (1 - c[0]) * (2 * c[1] + c[0]);
  dnds[9]  = .25 * (1 + c[0]) * (2 * c[1] - c[0]);
  dnds[10] = .25 * (1 + c[0]) * (2 * c[1] + c[0]);
  dnds[11] = .25 * (1 - c[0]) * (2 * c[1] - c[0]);
  dnds[12] = -   .5 * (1 - c[0] * c[0]);
  dnds[13] = - c[1] * (1 + c[0]);
  dnds[14] =     .5 * (1 - c[0] * c[0]);
  dnds[15] = - c[1] * (1 - c[0]);
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
template<> inline Real ElementClass<_quadrangle_4>::getInradius(const Real * coord) {
  Real a = Math::distance_2d(coord + 0, coord + 2);
  Real b = Math::distance_2d(coord + 2, coord + 4);
  Real c = Math::distance_2d(coord + 4, coord + 6);
  Real d = Math::distance_2d(coord + 6, coord + 0);

  Real h = std::min(a, std::min(b, std::min(c, d)));

  return h;
}
