/**
 * @file   element_class_point_1_inline_impl.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Feb 28 2020
 *
 * @brief  Specialization of the element_class class for the type _point_1
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/**
 * @verbatim
      x
    (0)
 @endverbatim
 *
 * @f{eqnarray*}{
 * N1 &=& 1
 * @f}
 *
 * @f{eqnarray*}{
 * q_0 &=& 0
 * @f}
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_point_1, _gt_point, _itp_lagrange_point_1,
                                     _ek_regular, 0, _git_point, 1);

/* --------------r------------------------------------------------------------
 */
template <>
inline void InterpolationElement<_itp_lagrange_point_1>::computeShapes(
    const Ref<const VectorXr> & /*natural_coords*/, Ref<VectorXr> N) {
  N(0) = 1; /// N1(q_0)
}
/* -------------------------------------------------------------------------- */
template <>
inline void InterpolationElement<_itp_lagrange_point_1>::computeDNDS(
    const Ref<const VectorXr> & /*natural_coords*/, Ref<MatrixXr> /*dnds*/) {}

/* -------------------------------------------------------------------------- */
template <>
inline void InterpolationElement<_itp_lagrange_point_1>::computeSpecialJacobian(
    const Ref<const MatrixXr> & /*J*/, Real & jac) {
  jac = 0.;
}

/* -------------------------------------------------------------------------- */
template <>
inline Real GeometricalElement<_gt_point>::getInradius(
    const Ref<const MatrixXr> & /*coord*/) {
  return 0.;
}
} // namespace akantu
