/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 * @section DESCRIPTION
 *
 * @verbatim
       \eta
     ^
     |
     x (0,0,1)
     |`
     |  `
     |` q `
     |  ` ° `
     x--------x----->  \xi
    (1,0,0)      (0,1,0)
 @endverbatim
 *
 * @subsection shapes shape functions
 * Parent:
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta \\
 * N2 &=& \xi \\
 * N3 &=& \eta
 * @f}
 * Sub 1:
 * @f{eqnarray*}{
 * N2 &=& \xi \\
 * N3 &=& \eta
 * @f}
 * Sub 2:
 * @f[
 * \begin{array}{lll}
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
 * \xi_{q0}  &=& 1/3 \qquad  \eta_{q0} = 1/3
 * @f}
 */

/* -------------------------------------------------------------------------- */
#include "element_class_igfem.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_IGFEM_ELEMENT_CLASS_PROPERTY(_igfem_triangle_4,
                                           _gt_igfem_triangle_4,
                                           _itp_igfem_triangle_4, _triangle_3,
                                           _triangle_3, _triangle_3, _ek_igfem,
                                           2, 1);

/* -------------------------------------------------------------------------- */
template <>
inline UInt ElementClass<_igfem_triangle_4>::getOrientation(
    const Vector<bool> & is_inside) {
  UInt orientation = 0;
  if (is_inside(1) && !is_inside(2))
    orientation = 0;

  else if (!is_inside(1) && is_inside(2))
    orientation = 1;

  else if (is_inside(1) && is_inside(2))
    orientation = 2;

  else if (!is_inside(1) && !is_inside(2))
    orientation = 3;

  return orientation;
}
} // namespace akantu
