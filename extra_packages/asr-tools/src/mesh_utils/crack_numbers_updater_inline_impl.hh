/**
 * @file   crack_numbers_updater_inline_impl.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date Tue Feb 8  2022
 * @brief  inline implementation for the crack number synchronizer
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
#include "communicator.hh"
#include "crack_numbers_updater.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CRACK_NUMBERS_UPDATER_INLINE_IMPL_CC__
#define __AKANTU_CRACK_NUMBERS_UPDATER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt
CrackNumbersUpdater::getNbData(const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  UInt size = 0;
  if (tag == SynchronizationTag::_crack_nb) {
    size += sizeof(int);
    auto & mesh = model.getMesh();
    for (auto elements_range : MeshElementsByTypes(elements)) {
      auto type = elements_range.getType();
      if (mesh.getKind(type) == _ek_cohesive) {
        UInt nb_elements = elements_range.getElements().size();
        size += nb_elements * sizeof(UInt);
      }
    }
  }
  return size;
}

/* -------------------------------------------------------------------------- */
inline void
CrackNumbersUpdater::packData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              const SynchronizationTag & tag) const {
  if (tag != SynchronizationTag::_crack_nb) {
    return;
  }

  auto & mesh = this->model.getMesh();
  auto && comm = mesh.getCommunicator();
  int prank = comm.whoAmI();
  buffer << prank;

  for (auto elements_range : MeshElementsByTypes(elements)) {
    auto type = elements_range.getType();
    auto gt = elements_range.getGhostType();
    if ((mesh.getKind(type) == _ek_cohesive) and (gt == _not_ghost)) {
      auto & coh_el_ids = elements_range.getElements();

      for (UInt coh_el : coh_el_ids) {
        auto & crack_numbers =
            model.getMesh().getData<UInt>("crack_numbers", type);
        buffer << crack_numbers(coh_el);
      }
    }
  }
}
/* --------------------------------------------------------------------------
 */
inline void CrackNumbersUpdater::unpackData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            const SynchronizationTag & tag) {
  if (tag != SynchronizationTag::_crack_nb) {
    return;
  }

  int proc;
  buffer >> proc;

  auto & mesh = model.getMesh();
  for (auto elements_range : MeshElementsByTypes(elements)) {
    auto type = elements_range.getType();
    auto gt = elements_range.getGhostType();
    if ((mesh.getKind(type) == _ek_cohesive) and (gt == _ghost)) {
      if (not mesh.hasData<UInt>("crack_numbers", type, gt)) {
        mesh.getDataPointer<UInt>("crack_numbers", type, gt);
      }

      auto & crack_numbers = mesh.getData<UInt>("crack_numbers", type, gt);

      auto & ghost_coh_el_ids = elements_range.getElements();
      for (UInt ghost_coh_el : ghost_coh_el_ids) {
        // resize crack_nb array on each iteration
        if (crack_numbers.size() <= ghost_coh_el) {
          crack_numbers.resize(ghost_coh_el + 1);
        }

        buffer >> crack_numbers(ghost_coh_el);
      }
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_CRACK_NUMBERS_UPDATER_INLINE_IMPL_CC__ */
