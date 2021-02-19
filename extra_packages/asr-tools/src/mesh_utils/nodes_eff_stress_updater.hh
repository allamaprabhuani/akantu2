/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_NODES_EFF_STRESS_UPDATER_HH__
#define __AKANTU_NODES_EFF_STRESS_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
class ElementSynchronizer;
}

namespace akantu {

class NodesEffStressUpdater : public DataAccessor<Element> {
public:
  NodesEffStressUpdater(Mesh & mesh, ElementSynchronizer & synchronizer,
                        Array<Real> & nodes_eff_stress)
      : mesh(mesh), synchronizer(synchronizer),
        nodes_eff_stress(nodes_eff_stress) {
    AKANTU_DEBUG_ASSERT(mesh.getNbNodes() == nodes_eff_stress.size(),
                        "Array nodes_eff_stress does not have the same size as "
                        "the number of nodes in the mesh");
  }

  void updateMaxEffStressAtNodes();

  /* ----------------------------------------------------------------- */
  /* Data Accessor inherited members                                   */
  /* ----------------------------------------------------------------- */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ----------------------------------------------------------------- */
  /* Members                                                           */
  /* ----------------------------------------------------------------- */
private:
  /// Reference to the mesh to update
  Mesh & mesh;

  /// distributed synchronizer to communicate nodes positions
  ElementSynchronizer & synchronizer;

  /// average effective stress of facets loop around a node
  Array<Real> & nodes_eff_stress;
};

} // namespace akantu

#include "nodes_eff_stress_updater_inline_impl.hh"

#endif /* __AKANTU_NODES_EFF_STRESS_UPDATER_HH__ */
