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
  if (tag == SynchronizationTag::_asr) {
    size += sizeof(int);
    // size += 2 * sizeof(Real);
    auto & mesh = model.getMesh();
    for (auto elements_range : MeshElementsByTypes(elements)) {
      auto type = elements_range.getType();
      if (mesh.getKind(type) == _ek_cohesive) {
        UInt nb_elements = elements_range.getElements().size();
        UInt nb_quad_cohesive =
            model.getFEEngine("CohesiveFEEngine").getNbIntegrationPoints(type);
        size += nb_elements * nb_quad_cohesive * sizeof(Real);
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
  if (tag != SynchronizationTag::_asr)
    return;

  int prank = model.getMesh().getCommunicator().whoAmI();
  buffer << prank;

  auto & mesh = model.getMesh();
  for (auto elements_range : MeshElementsByTypes(elements)) {
    auto type = elements_range.getType();
    auto gt = elements_range.getGhostType();
    if ((mesh.getKind(type) == _ek_cohesive) and (gt == _not_ghost)) {
      UInt nb_quad_cohesive =
          model.getFEEngine("CohesiveFEEngine").getNbIntegrationPoints(type);
      auto & coh_el_ids = elements_range.getElements();

      for (UInt coh_el : coh_el_ids) {
        const Array<UInt> & material_index_vec =
            model.getMaterialByElement(type, gt);
        Material & material = model.getMaterial(material_index_vec(coh_el));
        auto & crack_nbs = material.getInternal<Real>("crack_number")(type, gt);
        auto crack_nbs_it = crack_nbs.begin();
        for (auto i : arange(nb_quad_cohesive)) {
          buffer << crack_nbs_it[coh_el * nb_quad_cohesive + i];
        }
      }
    }
  }
}
/* --------------------------------------------------------------------------
 */
inline void CrackNumbersUpdater::unpackData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            const SynchronizationTag & tag) {
  if (tag != SynchronizationTag::_asr)
    return;

  int proc;
  buffer >> proc;

  auto & mesh = model.getMesh();
  for (auto elements_range : MeshElementsByTypes(elements)) {
    auto type = elements_range.getType();
    auto gt = elements_range.getGhostType();
    if ((mesh.getKind(type) == _ek_cohesive) and (gt == _ghost)) {
      UInt nb_quad_cohesive =
          model.getFEEngine("CohesiveFEEngine").getNbIntegrationPoints(type);
      auto & ghost_coh_el_ids = elements_range.getElements();
      for (UInt ghost_coh_el : ghost_coh_el_ids) {
        const Array<UInt> & material_index_vec =
            model.getMaterialByElement(type, gt);
        Material & material =
            model.getMaterial(material_index_vec(ghost_coh_el));
        auto & crack_nbs = material.getInternal<Real>("crack_number")(type, gt);
        auto crack_nbs_it = crack_nbs.begin();
        for (auto i : arange(nb_quad_cohesive)) {
          Real unpacked_crack_nb;
          buffer >> unpacked_crack_nb;
          crack_nbs_it[ghost_coh_el * nb_quad_cohesive + i] = unpacked_crack_nb;
        }
      }
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_CRACK_NUMBERS_UPDATER_INLINE_IMPL_CC__ */
