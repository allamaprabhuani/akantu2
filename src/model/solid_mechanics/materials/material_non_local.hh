/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 */

/* -------------------------------------------------------------------------- */
#include "constitutive_law_non_local_interface.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_NON_LOCAL_HH_
#define AKANTU_MATERIAL_NON_LOCAL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
class MaterialNonLocalInterface : public ConstitutiveLawNonLocalInterface {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law
  virtual void computeNonLocalStresses(GhostType ghost_type = _not_ghost) = 0;
};
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <Int dim, class LocalParent>
using MaterialNonLocal =
    ConstitutiveLawNonLocal<dim, MaterialNonLocalInterface, LocalParent>;

} // namespace akantu

#endif /* AKANTU_MATERIAL_NON_LOCAL_HH_ */
