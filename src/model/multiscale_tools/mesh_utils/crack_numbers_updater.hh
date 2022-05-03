/**
 * @file   crack_numbers_updater.hh
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

#ifndef __AKANTU_CRACK_NUMBERS_UPDATER_HH__
#define __AKANTU_CRACK_NUMBERS_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
class ElementSynchronizer;
} // namespace akantu

namespace akantu {

class CrackNumbersUpdater : public DataAccessor<Element> {
public:
  CrackNumbersUpdater(SolidMechanicsModel & model,
                      const ElementSynchronizer & synchronizer)
      : model(model), synchronizer(synchronizer) {}

  void communicateCrackNumbers();

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// Reference to the model
  SolidMechanicsModel & model;

  /// distributed synchronizer to communicate nodes positions
  const ElementSynchronizer & synchronizer;
};

} // namespace akantu

#include "crack_numbers_updater_inline_impl.hh"

#endif /* __AKANTU_CRACK_NUMBERS_UPDATER_HH__ */
