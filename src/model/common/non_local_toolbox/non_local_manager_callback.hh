/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH_
#define AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH_

namespace akantu {
class NonLocalManager;
} // namespace akantu

namespace akantu {

class NonLocalManagerCallback {
public:
  NonLocalManagerCallback() = default;
  NonLocalManagerCallback(const NonLocalManagerCallback &) = default;
  NonLocalManagerCallback(NonLocalManagerCallback &&) = default;
  NonLocalManagerCallback &
  operator=(const NonLocalManagerCallback &) = default;
  NonLocalManagerCallback & operator=(NonLocalManagerCallback &&) = default;

  virtual ~NonLocalManagerCallback() = default;

  virtual void initializeNonLocal() {}

  /* ------------------------------------------------------------------------ */
  virtual void
  insertIntegrationPointsInNeighborhoods(GhostType /*ghost_type*/) {
    AKANTU_TO_IMPLEMENT();
  }

  virtual void computeNonLocalContribution(GhostType /*ghost_type*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// update the values of the non local internal
  virtual void updateLocalInternal(ElementTypeMapReal & /*internal_flat*/,
                                   GhostType /*ghost_type*/,
                                   ElementKind /*kind*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// copy the results of the averaging in the materials
  virtual void updateNonLocalInternal(ElementTypeMapReal & /*internal_flat*/,
                                      GhostType /*ghost_type*/,
                                      ElementKind /*kind*/) {
    AKANTU_TO_IMPLEMENT();
  }
};

} // namespace akantu

#endif /* AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH_ */
