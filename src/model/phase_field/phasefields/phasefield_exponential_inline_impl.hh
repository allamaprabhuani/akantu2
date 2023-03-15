/**
 * Copyright (©) 2022-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 * 
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 * 
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

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

  for (Int i = 0; i < spatial_dimension; i++) {
    strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
    strain_diag_minus(i, i) = std::min(Real(0.), strain_values(i));
  }

  Matrix<Real> mat_tmp(spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_minus(spatial_dimension, spatial_dimension);

  strain_plus = strain_dir * strain_diag_plus * strain_dir.transpose();
  strain_minus = strain_dir * strain_diag_minus * strain_dir.transpose();

  trace_plus = std::max(Real(0.), strain_quad.trace());
  trace_minus = std::min(Real(0.), strain_quad.trace());

  for (Int i = 0; i < spatial_dimension; i++) {
    for (Int j = 0; j < spatial_dimension; j++) {
      sigma_plus(i, j) = static_cast<Real>(i == j) * lambda * trace_plus +
                         2 * mu * strain_plus(i, j);
      sigma_minus(i, j) = static_cast<Real>(i == j) * lambda * trace_minus +
                          2 * mu * strain_minus(i, j);
    }
  }

  phi_quad = sigma_plus.doubleDot(strain_plus) / 2.;
  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

} // namespace akantu
