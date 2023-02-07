#include "phasefield_exponential.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void
PhaseFieldExponential::computeDrivingForceOnQuad(const Real & phi_quad,
                                                 Real & driving_force_quad) {
  driving_force_quad = 2.0 * phi_quad;
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

  Matrix<Real> strain_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_minus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_dir(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_minus(spatial_dimension, spatial_dimension);

  Vector<Real> strain_values(spatial_dimension);

  Real trace_plus, trace_minus;

  strain_plus.zero();
  strain_minus.zero();
  strain_dir.zero();
  strain_values.zero();
  strain_diag_plus.zero();
  strain_diag_minus.zero();

  strain_quad.eig(strain_values, strain_dir);

  for (UInt i = 0; i < spatial_dimension; i++) {
    strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
    strain_diag_minus(i, i) = std::min(Real(0.), strain_values(i));
  }

  Matrix<Real> mat_tmp(spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_minus(spatial_dimension, spatial_dimension);

  mat_tmp.mul<false, true>(strain_diag_plus, strain_dir);
  strain_plus.mul<false, false>(strain_dir, mat_tmp);
  mat_tmp.mul<false, true>(strain_diag_minus, strain_dir);
  strain_minus.mul<false, true>(strain_dir, mat_tmp);

  trace_plus = std::max(Real(0.), strain_quad.trace());
  trace_minus = std::min(Real(0.), strain_quad.trace());

  for (UInt i = 0; i < spatial_dimension; i++) {
    for (UInt j = 0; j < spatial_dimension; j++) {
      sigma_plus(i, j) = static_cast<double>(i == j) * lambda * trace_plus +
                         2 * mu * strain_plus(i, j);
      sigma_minus(i, j) = static_cast<double>(i == j) * lambda * trace_minus +
                          2 * mu * strain_minus(i, j);
    }
  }

  phi_quad = 0.5 * sigma_plus.doubleDot(strain_plus);
  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

} //namespace akantu
