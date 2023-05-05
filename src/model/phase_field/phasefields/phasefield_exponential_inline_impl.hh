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
#include <algorithm>

namespace akantu {

inline void PhaseFieldExponential::computeDissipatedEnergyOnQuad(
    const Real & dam, const Vector<Real> & grad_d, Real & edis,
    Real & g_c_quad) {

  edis = 0.;
  for (auto i : arange(spatial_dimension)) {
    edis += 0.5 * g_c_quad * this->l0 * grad_d[i] * grad_d[i];
  }

  edis += g_c_quad * dam * dam / (2. * this->l0);
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

  Real trace = strain_quad.trace();
  Real trace_plus = std::max(Real(0.), trace);

  Matrix<Real> strain_dev(spatial_dimension, spatial_dimension);
  strain_dev = strain_quad -
               trace / Real(spatial_dimension) *
                   Matrix<Real>::Identity(spatial_dimension, spatial_dimension);

  Real kpa = this->lambda + 2. * this->mu / Real(spatial_dimension);

  phi_quad = 0.5 * kpa * trace_plus * trace_plus +
             this->mu * strain_dev.doubleDot(strain_dev);

  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

/* -------------------------------------------------------------------------- */
inline void PhaseFieldExponential::computePhiIsotropicOnQuad(
    const Matrix<Real> & strain_quad, Real & phi_quad, Real & phi_hist_quad) {

  Real trace = strain_quad.trace();

  phi_quad = 0.5 * this->lambda * trace * trace +
             this->mu * strain_quad.doubleDot(strain_quad);

  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

} // namespace akantu
