/* ------------------------------------------------------------------ */
#include "nodes_eff_stress_updater.hh"
#include "element_synchronizer.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* ------------------------------------------------------------------- */
#include <numeric>
/* ------------------------------------------------------------------- */

namespace akantu {

void NodesEffStressUpdater::updateMaxEffStressAtNodes() {
  if (mesh.getCommunicator().getNbProc() == 1)
    return;
  this->synchronizer.synchronizeOnce(*this, SynchronizationTag::_border_nodes);
  this->synchronizer.slaveReductionOnce(*this,
                                        SynchronizationTag::_border_nodes);
}

} // namespace akantu
