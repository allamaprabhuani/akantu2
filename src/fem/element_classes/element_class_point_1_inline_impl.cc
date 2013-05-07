/**
 * @file   element_class_point_1_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri May 2 09:09:21 2013
 *
 * @brief  Specialization of the element_class class for the type _point_1
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
	  x
    (0)
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * N1 &=& 1
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \q0 &=& 0
 * @f}
 */

AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_point_1, _gt_point, _itp_lagrange_point_1, _ek_regular, 0);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void
InterpolationElement<_itp_lagrange_point_1>::computeShapes(__attribute__ ((unused)) const vector_type & natural_coords,
                                                             vector_type & N) {
  N(0) = 1; /// N1(q_0)
}
/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void
InterpolationElement<_itp_lagrange_point_1>::computeDNDS(__attribute__ ((unused)) const vector_type & natural_coords,
                                                         __attribute__ ((unused)) matrix_type & dnds) {}


/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_lagrange_point_1>::computeSpecialJacobian(__attribute__ ((unused)) const Matrix<Real> & J,
                                                                    Real & jac){
   jac = 0.;
}



/* -------------------------------------------------------------------------- */
template<>
inline Real
GeometricalElement<_gt_point>::getInradius(const Matrix<Real> & coord) {
  return 0.;
}
