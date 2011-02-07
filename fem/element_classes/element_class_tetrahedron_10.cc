/**
 * @file   element_class_tetrahedron_10.cc
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @date   Thu Dec 1 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _tetrahedron_10
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
	\zeta
	  ^
	  |
       (0,0,1)
	  x
	  |` .
	  |  `  .
	  |    `   .
	  |      `    .  (0,0.5,0.5)
	  |        `    x.
	  |     q4 o `      .                   \eta
	  |            `       .             -,
(0,0,0.5) x             ` x (0.5,0,0.5)  -
	  |                `        x-(0,1,0)
	  |              q3 o`   -   '
	  |       (0,0.5,0)  - `      '
	  |             x-       `     x (0.5,0.5,0)
	  |     q1 o -         o q2`    '
	  |      -                   `   '
	  |  -                         `  '
	  x---------------x--------------` x-----> \xi
       (0,0,0)        (0.5,0,0)        (1,0,0)
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{lll}
 *   \xi_{0}  = 0   &  \eta_{0}  = 0   &  \zeta_{0}  = 0   \\
 *   \xi_{1}  = 1   &  \eta_{1}  = 0   &  \zeta_{1}  = 0   \\
 *   \xi_{2}  = 0   &  \eta_{2}  = 1   &  \zeta_{2}  = 0   \\
 *   \xi_{3}  = 0   &  \eta_{3}  = 0   &  \zeta_{3}  = 1   \\
 *   \xi_{4}  = 1/2 &  \eta_{4}  = 0   &  \zeta_{4}  = 0   \\
 *   \xi_{5}  = 1/2 &  \eta_{5}  = 1/2 &  \zeta_{5}  = 0   \\
 *   \xi_{6}  = 0   &  \eta_{6}  = 1/2 &  \zeta_{6}  = 0   \\
 *   \xi_{7}  = 0   &  \eta_{7}  = 0   &  \zeta_{7}  = 1/2 \\
 *   \xi_{8}  = 1/2 &  \eta_{8}  = 0   &  \zeta_{8}  = 1/2 \\
 *   \xi_{9}  = 0   &  \eta_{9}  = 1/2 &  \zeta_{9}  = 1/2
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{llll}
 *     N1  = (1 - \xi - \eta - \zeta) (1 - 2 \xi - 2 \eta - 2 \zeta)
 *           & \frac{\partial N1}{\partial \xi}    = 4 \xi + 4 \eta + 4 \zeta - 3
 *           & \frac{\partial N1}{\partial \eta}   = 4 \xi + 4 \eta + 4 \zeta - 3
 *           & \frac{\partial N1}{\partial \zeta}  = 4 \xi + 4 \eta + 4 \zeta - 3 \\
 *     N2  = \xi (2 \xi - 1)
 *           & \frac{\partial N2}{\partial \xi}    = 4 \xi - 1
 *           & \frac{\partial N2}{\partial \eta}   = 0
 *           & \frac{\partial N2}{\partial \zeta}  = 0 \\
 *     N3  = \eta (2 \eta - 1)
 *           & \frac{\partial N3}{\partial \xi}    = 0
 *           & \frac{\partial N3}{\partial \eta}   = 4 \eta - 1
 *           & \frac{\partial N3}{\partial \zeta}  = 0 \\
 *     N4  = \zeta (2 \zeta - 1)
 *           & \frac{\partial N4}{\partial \xi}    = 0
 *           & \frac{\partial N4}{\partial \eta}   = 0
 *           & \frac{\partial N4}{\partial \zeta}  = 4 \zeta - 1 \\
 *     N5  = 4 \xi (1 - \xi - \eta - \zeta)
 *           & \frac{\partial N5}{\partial \xi}    = 4 - 8 \xi - 4 \eta - 4 \zeta
 *           & \frac{\partial N5}{\partial \eta}   = -4 \xi
 *           & \frac{\partial N5}{\partial \zeta}  = -4 \xi \\
 *     N6  = 4 \xi \eta
 *           & \frac{\partial N6}{\partial \xi}    = 4 \eta
 *           & \frac{\partial N6}{\partial \eta}   = 4 \xi
 *           & \frac{\partial N6}{\partial \zeta}  = 0 \\
 *     N7  = 4 \eta (1 - \xi - \eta - \zeta)
 *           & \frac{\partial N7}{\partial \xi}    = -4 \eta
 *           & \frac{\partial N7}{\partial \eta}   = 4 - 4 \xi - 8 \eta - 4 \zeta
 *           & \frac{\partial N7}{\partial \zeta}  = -4 \eta \\
 *     N8  = 4 \zeta (1 - \xi - \eta - \zeta)
 *           & \frac{\partial N8}{\partial \xi}    = -4 \zeta
 *           & \frac{\partial N8}{\partial \eta}   = -4 \zeta
 *           & \frac{\partial N8}{\partial \zeta}  = 4 - 4 \xi - 4 \eta - 8 \zeta \\
 *     N9  = 4 \zeta \xi
 *           & \frac{\partial N9}{\partial \xi}    = 4 \zeta
 *           & \frac{\partial N9}{\partial \eta}   = 0
 *           & \frac{\partial N9}{\partial \zeta}  = 4 \xi \\
 *     N10 = 4 \eta \zeta
 *           & \frac{\partial N10}{\partial \xi}   = 0
 *           & \frac{\partial N10}{\partial \eta}  = 4 \zeta
 *           & \frac{\partial N10}{\partial \zeta} = 4 \eta \\
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f[
 * a = \frac{5 - \sqrt{5}}{20}\\
 * b = \frac{5 + 3 \sqrt{5}}{20}
 * \begin{array}{lll}
 *   \xi_{q_0}  = a   &  \eta_{q_0}  = a   &  \zeta_{q_0}  = a \\
 *   \xi_{q_1}  = b   &  \eta_{q_1}  = a   &  \zeta_{q_1}  = a \\
 *   \xi_{q_2}  = a   &  \eta_{q_2}  = b   &  \zeta_{q_2}  = a \\
 *   \xi_{q_3}  = a   &  \eta_{q_3}  = a   &  \zeta_{q_3}  = b
 * \end{array}
 * @f]
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedron_10>::nb_nodes_per_element;
template<> UInt ElementClass<_tetrahedron_10>::nb_quadrature_points;
template<> UInt ElementClass<_tetrahedron_10>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_tetrahedron_10>::computeShapes(const Real * natural_coords, 
								   Real * shapes){
  /// Natural coordinates
  Real xi = natural_coords[0];
  Real eta = natural_coords[1];
  Real zeta = natural_coords[2];
  Real sum = xi + eta + zeta;
  Real c0  = 1 - sum; 
  Real c1  = 1 - 2*sum;
  Real c2  = 2*xi - 1;
  Real c3  = 2*eta - 1;
  Real c4  = 2*zeta - 1;

  /// Shape functions
  shapes[0] = c0  * c1;
  shapes[1] = xi  * c2;
  shapes[2] = eta  * c3;
  shapes[3] = zeta * c4;
  shapes[4] = 4 * xi  * c0;
  shapes[5] = 4 * xi * eta;
  shapes[6] = 4 * eta * c0;
  shapes[7] = 4 * zeta * c0;
  shapes[8] = 4 * xi * zeta;
  shapes[9] = 4 * eta * zeta;
}

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_tetrahedron_10>::computeDNDS(__attribute__ ((unused)) const Real * natural_coords,
								 Real * dnds) {

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccccccc}
   *            \frac{\partial N1}{\partial \xi}   & \frac{\partial N2}{\partial \xi}
   *          & \frac{\partial N3}{\partial \xi}   & \frac{\partial N4}{\partial \xi}
   *          & \frac{\partial N5}{\partial \xi}   & \frac{\partial N6}{\partial \xi}
   *          & \frac{\partial N7}{\partial \xi}   & \frac{\partial N8}{\partial \xi}
   *          & \frac{\partial N9}{\partial \xi}   & \frac{\partial N10}{\partial \xi} \\
   *            \frac{\partial N1}{\partial \eta}  & \frac{\partial N2}{\partial \eta}
   *          & \frac{\partial N3}{\partial \eta}  & \frac{\partial N4}{\partial \eta}
   *          & \frac{\partial N5}{\partial \eta}  & \frac{\partial N6}{\partial \eta}
   *          & \frac{\partial N7}{\partial \eta}  & \frac{\partial N8}{\partial \eta}
   *          & \frac{\partial N9}{\partial \eta}  & \frac{\partial N10}{\partial \eta} \\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial \zeta}
   *          & \frac{\partial N3}{\partial \zeta} & \frac{\partial N4}{\partial \zeta}
   *          & \frac{\partial N5}{\partial \zeta} & \frac{\partial N6}{\partial \zeta}
   *          & \frac{\partial N7}{\partial \zeta} & \frac{\partial N8}{\partial \zeta}
   *          & \frac{\partial N9}{\partial \zeta} & \frac{\partial N10}{\partial \zeta}
   *          \end{array}
   *        \right)
   * @f]
   */

  /// Natural coordinates
  Real xi = natural_coords[0];
  Real eta = natural_coords[1];
  Real zeta = natural_coords[2];
  Real sum = xi + eta + zeta;

  /// dN/dxi
  dnds[0] = 4 * sum - 3;
  dnds[1] = 4 * xi - 1;
  dnds[2] = 0;
  dnds[3] = 0;
  dnds[4] = 4 * (1 - sum - xi);
  dnds[5] = 4 * eta;
  dnds[6] = -4 * eta;
  dnds[7] = -4 * zeta;
  dnds[8] = 4 * zeta;
  dnds[9] = 0;

  /// dN/deta
  dnds[10] = 4 * sum - 3;
  dnds[11] = 0;
  dnds[12] = 4 * eta - 1;
  dnds[13] = 0;
  dnds[14] = -4 * xi;
  dnds[15] = 4 * xi;
  dnds[16] = 4 * (1 - sum - eta);
  dnds[17] = -4 * zeta;
  dnds[18] = 0;
  dnds[19] = 4 * zeta;

  /// dN/dzeta
  dnds[20] = 4 * sum - 3;
  dnds[21] = 0;
  dnds[22] = 0;
  dnds[23] = 4 * zeta - 1;
  dnds[24] = -4 * xi;
  dnds[25] = 0;
  dnds[26] = -4 * eta;
  dnds[27] = 4 * (1 - sum - zeta);
  dnds[28] = 4 * xi;
  dnds[29] = 4 * eta;

}

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_tetrahedron_10>::computeJacobian(const Real * dxds,
								     const UInt dimension, 
								     Real & jac) {

  if (dimension == spatial_dimension){
    Real weight = 1./24.; 
    Real det_dxds = Math::det3(dxds);
    jac = det_dxds * weight;    
  }  
  else {
    AKANTU_DEBUG_ERROR("to be implemented");
  }
}
 
/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_tetrahedron_10>::getInradius(const Real * coord) {

  // Only take the four corner tetrahedra
  UInt tetrahedra[4][4] = {
    {0, 4, 6, 7},
    {4, 1, 5, 8},
    {6, 5, 2, 9},
    {7, 8, 9, 3}
  };

  Real inradius = std::numeric_limits<Real>::max();
  for (UInt t = 0; t < 4; t++) {
    Real ir = Math::tetrahedron_inradius(coord + tetrahedra[t][0] * spatial_dimension,
					 coord + tetrahedra[t][1] * spatial_dimension,
					 coord + tetrahedra[t][2] * spatial_dimension,
					 coord + tetrahedra[t][3] * spatial_dimension);
    inradius = ir < inradius ? ir : inradius;
  }

  return inradius;
}
