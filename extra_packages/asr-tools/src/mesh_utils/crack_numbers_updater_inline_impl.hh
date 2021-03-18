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
      // // resize crack_nb array if size of ghost coh elements is >
      // if (crack_numbers.size() != ghost_coh_el_ids.size())
      //   crack_numbers.resize(ghost_coh_el_ids.size());

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
