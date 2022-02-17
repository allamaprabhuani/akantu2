/**
 * @file   resolution_utils.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon May 20 2019
 * @date last modification: Sun Jun 06 2021
 *
 * @brief  All resolution utils necessary for various tasks
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "contact_element.hh"
#include "contact_mechanics_model.hh"
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
  /// computes the shape function matric for the contact element (@f$A
  /// @f$) where row is equal to spatial dimension and cols is equal
  /// to spatial dimension times number of nodes in contact element
  static void computeShapeFunctionMatric(const ContactElement & /*element*/,
                                         const Vector<Real> & /*projection*/,
                                         Matrix<Real> & /*shape_matric*/);
};

} // namespace akantu

#endif /* __AKANTU_RESOLUTION_UTILS_HH__ */
