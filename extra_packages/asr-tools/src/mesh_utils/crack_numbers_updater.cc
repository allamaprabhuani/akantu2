/* -------------------------------------------------------------------------- */
#include "crack_numbers_updater.hh"
#include "element_synchronizer.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

void CrackNumbersUpdater::communicateCrackNumbers() {
  if (model.getMesh().getCommunicator().getNbProc() == 1)
    return;
  this->synchronizer.synchronizeOnce(*this, SynchronizationTag::_asr);
}

} // namespace akantu
