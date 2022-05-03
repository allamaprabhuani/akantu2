/**
 * @file   crack_numbers_updater.cc
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date Tue Feb 8  2022
 * @brief  synchronizer for the crack numbers
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
#include "crack_numbers_updater.hh"
#include "element_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

void CrackNumbersUpdater::communicateCrackNumbers() {
  if (model.getMesh().getCommunicator().getNbProc() == 1)
    return;
  this->synchronizer.synchronizeOnce(*this, SynchronizationTag::_rve_crack_nb);
}

} // namespace akantu
