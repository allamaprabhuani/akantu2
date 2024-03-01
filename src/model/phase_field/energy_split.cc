
/* -------------------------------------------------------------------------- */
#include "energy_split.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
EnergySplit::EnergySplit(Real E, Real nu, bool plane_stress)
    : E(E), nu(nu), plane_stress(plane_stress) {
  this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  if (this->plane_stress) {
    this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - this->nu));
  }
  this->mu = this->E / (2 * (1 + this->nu));
}

} // namespace akantu
