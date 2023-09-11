#include "phasefield_linear.hh"

namespace akantu {

inline void
PhaseFieldLinear::computeDissipatedEnergyOnQuad(const Real & dam,
                                                const Vector<Real> & grad_d,
                                                Real & edis, Real & g_c_quad) {

  edis = 0.;
  for (auto i : arange(spatial_dimension)) {
    edis += 3. / 8 * g_c_quad * this->l0 * grad_d[i] * grad_d[i];
  }

  edis += 3 * g_c_quad * dam * dam / (8 * this->l0);
}

/* -------------------------------------------------------------------------- */
inline void PhaseFieldLinear::computeDrivingForceOnQuad(
    const Real & phi_quad, Real & driving_force_quad, const Real & g_c_quad) {
  driving_force_quad = phi_quad - 3 * g_c_quad / (16 * this->l0);
}

/* -------------------------------------------------------------------------- */
inline void PhaseFieldLinear::computeDamageEnergyDensityOnQuad(
    const Real & phi_quad, Real & dam_energy_quad, const Real & g_c_quad) {
  dam_energy_quad = 2 * phi_quad;
}

/* -------------------------------------------------------------------------- */
inline void PhaseFieldLinear::computePhiOnQuad(const Matrix<Real> & strain_quad,
                                               Real & phi_quad) {
  Real trace = strain_quad.trace();
  Real trace_plus = std::max(Real(0.), trace);

  Matrix<Real> strain_dev(spatial_dimension, spatial_dimension);
  strain_dev = strain_quad -
               trace / Real(spatial_dimension) *
                   Matrix<Real>::Identity(spatial_dimension, spatial_dimension);

  Real kpa = this->lambda + 2. * this->mu / Real(spatial_dimension);

  phi_quad = 0.5 * kpa * trace_plus * trace_plus +
             this->mu * strain_dev.doubleDot(strain_dev);
}

/* -------------------------------------------------------------------------- */
inline void
PhaseFieldLinear::computePhiIsotropicOnQuad(const Matrix<Real> & strain_quad,
                                            Real & phi_quad) {
  Real trace = strain_quad.trace();

  phi_quad = 0.5 * this->lambda * trace * trace +
             this->mu * strain_quad.doubleDot(strain_quad);
}

} // namespace akantu
