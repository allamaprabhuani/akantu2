/**
 * @file   element_class_hexahedron_8.cc
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @date   Mon Mar 14 17:27:44 2010
 *
 * @brief   Specialization of the element_class class for the type _hexahedron_8
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
                   \zeta
                    ^
         (-1,1,1)   |     (1,1,1)
                8---|------7
               /|   |     /|
              / |   |    / |
   (-1,-1,1) 5----------6  | (1,-1,1)
             |  |   |   |  |
             |  |   |   |  |
             |  |   +---|-------> \xi
             |  |  /    |  |
   (-1,1,-1) |  4-/-----|--3 (1,1,-1)
             | / /      | /
             |/ /       |/
             1-/--------2
   (-1,-1,-1) /        (1,-1,-1)
             /
            \eta
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{llll}
 * N1 = (1 - \xi) (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N1}{\partial \xi}  = - (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi) (1 - \zeta) / 8
 *       & \frac{\partial N1}{\partial \zeta} = - (1 - \xi) (1 - \eta) / 8 \\
 * N2 = (1 + \xi) (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta) (1 - \zeta) / 8
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi) (1 - \zeta) / 8
 *       & \frac{\partial N2}{\partial \zeta} = - (1 + \xi) (1 - \eta) / 8 \\
 * N3 = (1 + \xi) (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi) (1 - \zeta) / 8
 *       & \frac{\partial N3}{\partial \zeta} = - (1 + \xi) (1 + \eta) / 8 \\
 * N4 = (1 - \xi) (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N4}{\partial \xi}  = - (1 + \eta) (1 - \zeta) / 8
 *       & \frac{\partial N4}{\partial \eta} = (1 - \xi) (1 - \zeta) / 8
 *       & \frac{\partial N4}{\partial \zeta} = - (1 - \xi) (1 + \eta) / 8 \\
 * N5 = (1 - \xi) (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N5}{\partial \xi}  = - (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N5}{\partial \eta} = - (1 - \xi) (1 + \zeta) / 8
 *       & \frac{\partial N5}{\partial \zeta} = (1 - \xi) (1 - \eta) / 8 \\
 * N6 = (1 + \xi) (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N6}{\partial \xi}  = (1 - \eta) (1 + \zeta) / 8
 *       & \frac{\partial N6}{\partial \eta} = - (1 + \xi) (1 + \zeta) / 8
 *       & \frac{\partial N6}{\partial \zeta} = (1 + \xi) (1 - \eta) / 8 \\
 * N7 = (1 + \xi) (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N7}{\partial \xi}  = (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N7}{\partial \eta} = (1 + \xi) (1 + \zeta) / 8
 *       & \frac{\partial N7}{\partial \zeta} = (1 + \xi) (1 + \eta) / 8 \\
 * N8 = (1 - \xi) (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N8}{\partial \xi}  = - (1 + \eta) (1 + \zeta) / 8
 *       & \frac{\partial N8}{\partial \eta} = (1 - \xi) (1 + \zeta) / 8
 *       & \frac{\partial N8}{\partial \zeta} = (1 - \xi) (1 + \eta) / 8 \\
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& -1/\sqrt{3} \qquad  \eta_{q0} = -1/\sqrt{3} \qquad \zeta_{q0} = -1/\sqrt{3} \\
 * \xi_{q1}  &=&  1/\sqrt{3} \qquad  \eta_{q1} = -1/\sqrt{3} \qquad \zeta_{q1} = -1/\sqrt{3} \\
 * \xi_{q2}  &=&  1/\sqrt{3} \qquad  \eta_{q2} =  1/\sqrt{3} \qquad \zeta_{q2} = -1/\sqrt{3} \\
 * \xi_{q3}  &=& -1/\sqrt{3} \qquad  \eta_{q3} =  1/\sqrt{3} \qquad \zeta_{q3} = -1/\sqrt{3} \\
 * \xi_{q4}  &=& -1/\sqrt{3} \qquad  \eta_{q4} = -1/\sqrt{3} \qquad \zeta_{q4} =  1/\sqrt{3} \\
 * \xi_{q5}  &=&  1/\sqrt{3} \qquad  \eta_{q5} = -1/\sqrt{3} \qquad \zeta_{q5} =  1/\sqrt{3} \\
 * \xi_{q6}  &=&  1/\sqrt{3} \qquad  \eta_{q6} =  1/\sqrt{3} \qquad \zeta_{q6} =  1/\sqrt{3} \\
 * \xi_{q7}  &=& -1/\sqrt{3} \qquad  \eta_{q7} =  1/\sqrt{3} \qquad \zeta_{q7} =  1/\sqrt{3} \\
 * @f}
 */


/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_hexahedron_8>::nb_nodes_per_element;
template<> UInt ElementClass<_hexahedron_8>::nb_quadrature_points;
template<> UInt ElementClass<_hexahedron_8>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_hexahedron_8>::computeShapes(const Real * natural_coords,
								   Real * shapes) {
  /// Natural coordinates
  const Real * c = natural_coords;

  shapes[0] = .125 * (1 - c[0]) * (1 - c[1]) * (1 - c[3]); /// N1(q_0)
  shapes[1] = .125 * (1 + c[0]) * (1 - c[1]) * (1 - c[3]); /// N2(q_0)
  shapes[2] = .125 * (1 + c[0]) * (1 + c[1]) * (1 - c[3]); /// N3(q_0)
  shapes[3] = .125 * (1 - c[0]) * (1 + c[1]) * (1 - c[3]); /// N4(q_0)
  shapes[4] = .125 * (1 - c[0]) * (1 - c[1]) * (1 + c[3]); /// N5(q_0)
  shapes[5] = .125 * (1 + c[0]) * (1 - c[1]) * (1 + c[3]); /// N6(q_0)
  shapes[6] = .125 * (1 + c[0]) * (1 + c[1]) * (1 + c[3]); /// N7(q_0)
  shapes[7] = .125 * (1 - c[0]) * (1 + c[1]) * (1 + c[3]); /// N8(q_0)
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_hexahedron_8>::computeDNDS(const Real * natural_coords,
								 Real * dnds) {

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial \xi}
   *               & \frac{\partial N3}{\partial \xi}  & \frac{\partial N4}{\partial \xi}
   *               & \frac{\partial N5}{\partial \xi}  & \frac{\partial N6}{\partial \xi}
   *               & \frac{\partial N7}{\partial \xi}  & \frac{\partial N8}{\partial \xi}\\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial \eta}
   *               & \frac{\partial N3}{\partial \eta} & \frac{\partial N4}{\partial \eta}
   *               & \frac{\partial N5}{\partial \eta} & \frac{\partial N6}{\partial \eta}
   *               & \frac{\partial N7}{\partial \eta} & \frac{\partial N8}{\partial \eta}\\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial \zeta}
   *               & \frac{\partial N3}{\partial \zeta} & \frac{\partial N4}{\partial \zeta}
   *               & \frac{\partial N5}{\partial \zeta} & \frac{\partial N6}{\partial \zeta}
   *               & \frac{\partial N7}{\partial \zeta} & \frac{\partial N8}{\partial \zeta}
   *          \end{array}
   *        \right)
   * @f]
   */

  const Real * c = natural_coords;

  dnds[0]  = - .125 * (1 - c[1]) * (1 - c[2]);;
  dnds[1]  =   .125 * (1 - c[1]) * (1 - c[2]);;
  dnds[2]  =   .125 * (1 + c[1]) * (1 - c[2]);;
  dnds[3]  = - .125 * (1 + c[1]) * (1 - c[2]);;
  dnds[4]  = - .125 * (1 - c[1]) * (1 + c[2]);;
  dnds[5]  =   .125 * (1 - c[1]) * (1 + c[2]);;
  dnds[6]  =   .125 * (1 + c[1]) * (1 + c[2]);;
  dnds[7]  = - .125 * (1 + c[1]) * (1 + c[2]);;

  dnds[8]  = - .125 * (1 - c[0]) * (1 - c[2]);;
  dnds[9]  = - .125 * (1 + c[0]) * (1 - c[2]);;
  dnds[10] =   .125 * (1 + c[0]) * (1 - c[2]);;
  dnds[11] =   .125 * (1 - c[0]) * (1 - c[2]);;
  dnds[12] = - .125 * (1 - c[0]) * (1 + c[2]);;
  dnds[13] = - .125 * (1 + c[0]) * (1 + c[2]);;
  dnds[14] =   .125 * (1 + c[0]) * (1 + c[2]);;
  dnds[15] =   .125 * (1 - c[0]) * (1 + c[2]);;

  dnds[16] = - .125 * (1 - c[0]) * (1 - c[1]);;
  dnds[17] = - .125 * (1 + c[0]) * (1 - c[1]);;
  dnds[18] = - .125 * (1 + c[0]) * (1 + c[1]);;
  dnds[19] = - .125 * (1 - c[0]) * (1 + c[1]);;
  dnds[20] =   .125 * (1 - c[0]) * (1 - c[1]);;
  dnds[21] =   .125 * (1 + c[0]) * (1 - c[1]);;
  dnds[22] =   .125 * (1 + c[0]) * (1 + c[1]);;
  dnds[23] =   .125 * (1 - c[0]) * (1 + c[1]);;
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_hexahedron_8>::computeJacobian(const Real * dxds,
								   const UInt dimension,
								   Real & jac){
  Real weight = 1.; 
  if (dimension == spatial_dimension){
    Real det_dxds = Math::det3(dxds);
    jac = det_dxds * weight;
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}



/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_hexahedron_8>::getInradius(const Real * coord) {
  Real a = Math::distance_3d(coord +  0, coord +  3);
  Real b = Math::distance_3d(coord +  3, coord +  6);
  Real c = Math::distance_3d(coord +  6, coord +  9);
  Real d = Math::distance_3d(coord +  9, coord +  0);
  Real e = Math::distance_3d(coord +  0, coord + 12);
  Real f = Math::distance_3d(coord +  3, coord + 15);
  Real g = Math::distance_3d(coord +  6, coord + 18);
  Real h = Math::distance_3d(coord +  9, coord + 21);
  Real i = Math::distance_3d(coord + 12, coord + 15);
  Real j = Math::distance_3d(coord + 15, coord + 18);
  Real k = Math::distance_3d(coord + 18, coord + 21);
  Real l = Math::distance_3d(coord + 21, coord + 12);

  Real x = std::min(a, std::min(b, std::min(c, d)));
  Real y = std::min(e, std::min(f, std::min(g, h)));
  Real z = std::min(i, std::min(j, std::min(k, l)));
  Real p = std::min(x, std::min(y, z));

  return p;
}
