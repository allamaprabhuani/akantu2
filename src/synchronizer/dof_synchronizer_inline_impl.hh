/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communication_buffer.hh"
#include "data_accessor.hh"
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_CC__
// #define __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline Int DOFSynchronizer::canScatterSize() {
  return dof_manager.getLocalSystemSize();
}

/* -------------------------------------------------------------------------- */
inline Int DOFSynchronizer::gatheredSize() {
  return dof_manager.getSystemSize();
}

inline Idx DOFSynchronizer::localToGlobalEntity(const Idx & local) {
  return dof_manager.localToGlobalEquationNumber(local);
}
} // namespace akantu

//#endif /* __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_CC__ */
