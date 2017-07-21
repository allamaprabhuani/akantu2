/**
 * @file   non_local_manager_callback.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Sun Jul 09 2017
 *
 * @brief Callback functions for the non local manager
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
 */
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH__
#define __AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH__

namespace akantu {
class NonLocalManager;
} // namespace akantu

namespace akantu {

class NonLocalManagerCallback {
public:
  /* ------------------------------------------------------------------------ */
  virtual void
  insertIntegrationPointsInNeighborhoods(const GhostType & ghost_type) = 0;

  virtual void computeNonLocalStresses(const GhostType & ghost_type) = 0;

  /// update the values of the non local internal
  virtual void updateLocalInternal(ElementTypeMapReal & internal_flat,
                                   const GhostType & ghost_type,
                                   const ElementKind & kind) = 0;

  /// copy the results of the averaging in the materials
  virtual void updateNonLocalInternal(ElementTypeMapReal & internal_flat,
                                      const GhostType & ghost_type,
                                      const ElementKind & kind) = 0;
};

} // namespace akantu

#endif /* __AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH__ */
