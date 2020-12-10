/* -------------------------------------------------------------------------- */
#include "nodes_flag_updater.hh"
#include "element_synchronizer.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

void NodesFlagUpdater::fillPreventInsertion() {
  if (mesh.getCommunicator().getNbProc() == 1)
    return;
  this->synchronizer.slaveReductionOnce(*this,
                                        SynchronizationTag::_border_nodes);
}

} // namespace akantu
