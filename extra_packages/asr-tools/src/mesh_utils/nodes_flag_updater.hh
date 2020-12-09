/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_NODES_FLAG_UPDATER_HH__
#define __AKANTU_NODES_FLAG_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
class ElementSynchronizer;
} // namespace akantu

namespace akantu {

class NodesFlagUpdater : public DataAccessor<Element> {
public:
  NodesFlagUpdater(Mesh & mesh, ElementSynchronizer & synchronizer,
                   Array<bool> & prevent_insertion)
      : mesh(mesh), synchronizer(synchronizer),
        prevent_insertion(prevent_insertion) {
    AKANTU_DEBUG_ASSERT(
        mesh.getNbNodes() == prevent_insertion.size(),
        "Array prevent_insertion does not have the same size as "
        "the number of nodes in the mesh");
  }

  void fillPreventInsertion();

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
  /// Reference to the mesh to update
  Mesh & mesh;

  /// distributed synchronizer to communicate nodes positions
  ElementSynchronizer & synchronizer;

  /// Tells if a reduction is taking place or not
  bool reduce{false};

  /// tells at which nodes ASR segments shouldn't be inserted
  Array<bool> & prevent_insertion;
};

} // namespace akantu

#include "nodes_flag_updater_inline_impl.hh"

#endif /* __AKANTU_NODES_FLAG_UPDATER_HH__ */
