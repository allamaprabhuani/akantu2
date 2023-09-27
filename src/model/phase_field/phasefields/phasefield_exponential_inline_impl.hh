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
template <Int dim>
inline void PhaseFieldExponential<dim>::computeDissipatedEnergyOnQuad(
    const Real & dam, const Vector<Real> & grad_d, Real & edis,
    Real & g_c_quad) {

  edis = 0.;
  for (auto i : arange(spatial_dimension)) {
    edis += 0.5 * g_c_quad * this->l0 * grad_d[i] * grad_d[i];
  }

  edis += g_c_quad * dam * dam / (2 * this->l0);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldExponential<dim>::computeDamageEnergyDensityOnQuad(
    const Real & phi_quad, Real & dam_energy_quad, const Real & g_c_quad) {
  dam_energy_quad = 2.0 * phi_quad + g_c_quad / this->l0;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldExponential<dim>::computePhiOnQuad(
    const Matrix<Real> & strain_quad, Real & phi_quad, Real & phi_hist_quad) {
  Matrix<Real, dim, dim> strain_dir;
  Vector<Real, dim> strain_values;
  strain_quad.eig(strain_values, strain_dir);

  Matrix<Real, dim, dim> strain_diag_plus;
  strain_diag_plus.zero();
  for (Int i = 0; i < dim; i++) {
    strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
  }

  Matrix<Real, dim, dim> strain_plus =
      strain_dir * strain_diag_plus * strain_dir.transpose();

  auto trace_plus = std::max(0., strain_quad.trace());

  Matrix<Real, dim, dim> sigma_plus =
      Matrix<Real, dim, dim>::Identity() * lambda * trace_plus +
      2. * mu * strain_plus;

  phi_quad = sigma_plus.doubleDot(strain_plus) / 2.;
  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void PhaseFieldExponential<dim>::computePhiIsotropicOnQuad(
    const Matrix<Real> & strain_quad, Real & phi_quad, Real & phi_hist_quad) {

  Real trace = strain_quad.trace();

  phi_quad = 0.5 * this->lambda * trace * trace +
             this->mu * strain_quad.doubleDot(strain_quad);

  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

} // namespace akantu
