/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "nodes_flag_updater.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NODES_FLAG_UPDATER_INLINE_IMPL_CC__
#define __AKANTU_NODES_FLAG_UPDATER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt NodesFlagUpdater::getNbData(const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  UInt size = 0;
  if (tag == SynchronizationTag::_asr) {
    size +=
        Mesh::getNbNodesPerElementList(elements) * sizeof(bool) + sizeof(int);
  }
  return size;
}

/* -------------------------------------------------------------------------- */
inline void NodesFlagUpdater::packData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  if (tag != SynchronizationTag::_asr)
    return;

  int prank = mesh.getCommunicator().whoAmI();
  buffer << prank;

  for (auto & element : elements) {
    /// get element connectivity
    const Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      buffer << prevent_insertion(node);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void NodesFlagUpdater::unpackData(CommunicationBuffer & buffer,
                                         const Array<Element> & elements,
                                         const SynchronizationTag & tag) {
  if (tag != SynchronizationTag::_asr)
    return;

  MeshAccessor mesh_accessor(mesh);
  auto dim = mesh.getSpatialDimension();
  auto & nodes_pos = mesh_accessor.getNodes();
  auto pos_it = nodes_pos.begin(dim);

  int proc;
  buffer >> proc;

  for (auto & element : elements) {
    /// get element connectivity
    Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      bool unpacked_node_flag;
      buffer >> unpacked_node_flag;
      prevent_insertion(node) = true;
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_NODES_FLAG_UPDATER_INLINE_IMPL_CC__ */
