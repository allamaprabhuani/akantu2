/**
 * @file   resolution_utils.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon May 20 2019
 * @date last modification: Mon May 20 2019
 *
 * @brief  All resolution utils necessary for various tasks
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "contact_mechanics_model.hh"
#include "contact_element.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RESOLUTION_UTILS_HH__
#define __AKANTU_RESOLUTION_UTILS_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

class ResolutionUtils {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// computes the first covariant metric tensor (@f$A_{\alpha\beta}@f$) where @f$\alpha,
  /// \beta@f$ are surface directions
  static void computeMetricTensor(Matrix<Real> & metric_tensor,
                                  const Matrix<Real> & tangents);

  /// computes the second covariant metric tensor
  /// (@f$H_{\alpha\beta}@f$)
  static void computeSecondMetricTensor(const ContactElement &, const Matrix<Real> &,
					const Vector<Real> &, Matrix<Real> &);
  
  /// computes the first variation of normal gap
  static void firstVariationNormalGap(const ContactElement & , const Vector<Real> &,
				      const Vector<Real> &, Vector<Real> &);

  /// computes the seond variation of normal gap
  static void secondVariationNormalGap(const ContactElement & , const Matrix<Real> & ,
				       const Matrix<Real> &, const Vector<Real> &,
				       const Vector<Real> &, Real &,
				       Matrix<Real> & );
  
  /// computes (@f$N_{\alpha}@f$) where \alpha is surface dimension
  /// and it is shape derivatives times normal
  static void computeNalpha(const ContactElement & , const Vector<Real> &,
			    const Vector<Real> &, Array<Real> & );

  /// computes (@f$T_{\alpha}@f$) where @f$\alpha@f$ is surface
  /// dimension and it is shape functions times the tangents
  static void computeTalpha(const ContactElement &, const Matrix<Real> &,
			    const Vector<Real> &, Array<Real> & );

  /// computes (@f$\nabla \xi_{\alpha}@f$) where @f$\alpha@f$ is surface
  /// dimension
  static void firstVariationNaturalCoordinate(const ContactElement &, const Matrix<Real> &,
					      const Vector<Real> &, const Vector<Real> &,
					      const Real &, Array<Real> &);

  ///computes second variation of surface parameter
  static void secondvariationNaturalCoordinate(const Vector<Real> & projection,
					       const Vector<Real> & previous_projection,
					       const Element & element,
					       const Element & previous_element);
  
  /// computes @f$T_{\alpha\beta} @f$ which is shape derivatives
  /// times the tangents
  //static void computeTalphabeta(Array<Real> & t_alpha_beta,
  //                              ContactElement & element);

  /// computes @f$N_{\alpha\beta} @f$ which is shape 2nd derivatives times
  /// the normal
  //static void computeNalphabeta(Array<Real> & n_alpha_beta,
  //                              ContactElement & element);

  /// computes @f$P_{\alpha} @f$
  //static void computePalpha(Array<Real> & p_alpha, ContactElement & element);

  /// computes @f$G_{\alpha}@f$
  //static void computeGalpha(Array<Real> & g_alpha, Array<Real> & t_alpha_beta,
  //                         Array<Real> & d_alpha, Matrix<Real> & phi,
  //                          ContactElement &);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

#endif /* __AKANTU_RESOLUTION_UTILS_HH__ */
