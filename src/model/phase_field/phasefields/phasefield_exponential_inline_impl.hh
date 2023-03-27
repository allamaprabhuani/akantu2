#include "phasefield_exponential.hh"
#include <algorithm>

namespace akantu {

inline void PhaseFieldExponential::computeDissipatedEnergyOnQuad(
    const Real & dam, const Vector<Real> & grad_d, Real & edis,
    Real & g_c_quad) {

  edis = 0.;
  for (auto i : arange(spatial_dimension)) {
    edis += 0.5 * g_c_quad * this->l0 * grad_d[i] * grad_d[i];
  }

  edis += g_c_quad * dam * dam / (2 * this->l0);
}

/* -------------------------------------------------------------------------- */
inline void PhaseFieldExponential::computeDamageEnergyDensityOnQuad(
    const Real & phi_quad, Real & dam_energy_quad, const Real & g_c_quad) {
  dam_energy_quad = 2.0 * phi_quad + g_c_quad / this->l0;
}

/* -------------------------------------------------------------------------- */
inline void
PhaseFieldExponential::computePhiOnQuad(const Matrix<Real> & strain_quad,
                                        Real & phi_quad, Real & phi_hist_quad) {
  // Matrix<Real> strain_plus(spatial_dimension, spatial_dimension);
  // Matrix<Real> strain_dir(spatial_dimension, spatial_dimension);
  // Matrix<Real> strain_diag_plus(spatial_dimension, spatial_dimension);

  // Vector<Real> strain_values(spatial_dimension);

  Real trace = strain_quad.trace();
  Real trace_plus = std::max(Real(0.), trace);

  Matrix<Real> strain_dev(spatial_dimension, spatial_dimension);
  strain_dev = strain_quad -
               trace / static_cast<Real>(spatial_dimension) *
                   Matrix<Real>::Identity(spatial_dimension, spatial_dimension);

  Real Kn = this->lambda + 2. * this->mu / static_cast<Real>(spatial_dimension);

  phi_quad = 0.5 * Kn * trace_plus * trace_plus +
             this->mu * strain_dev.doubleDot(strain_dev);

  // strain_plus.zero();
  // strain_dir.zero();
  // strain_values.zero();
  // strain_diag_plus.zero();

  // strain_quad.eig(strain_values, strain_dir);

  // for (Int i = 0; i < spatial_dimension; i++) {
  //   strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
  // }

  // Matrix<Real> mat_tmp(spatial_dimension, spatial_dimension);
  // Matrix<Real> sigma_plus(spatial_dimension, spatial_dimension);

  // mat_tmp = strain_diag_plus * strain_dir.transpose();
  // strain_plus = strain_dir * mat_tmp;

  // trace_plus = std::max(Real(0.), strain_quad.trace());

  // for (Int i = 0; i < spatial_dimension; i++) {
  //   for (Int j = 0; j < spatial_dimension; j++) {
  //     sigma_plus(i, j) = static_cast<Real>(i == j) * lambda * trace_plus +
  //                        2 * mu * strain_plus(i, j);
  //   }
  // }

  // phi_quad = sigma_plus.doubleDot(strain_plus) / 2.;
  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

} // namespace akantu
