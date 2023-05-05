/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "material_phasefield.hh"
#include <algorithm>
#include <iostream>

#ifndef AKANTU_MATERIAL_PHASEFIELD_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_PHASEFIELD_INLINE_IMPL_HH_
/* -------------------------------------------------------------------------- */
namespace akantu {
template <Int dim>
template <class Args>
inline void MaterialPhaseField<dim>::computeStressOnQuad(Args && args) {
  MaterialElastic<dim>::computeStressOnQuad(args);

  auto && dam = args["damage"_n];
  args["sigma"_n] *= (1 - dam) * (1 - dam) + eta;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
void MaterialPhaseField<dim>::computeTangentModuliOnQuad(Args && args) {
  MaterialElastic<dim>::computeTangentModuliOnQuad(args);

  auto dam = args["damage"_n];
  args["tangent_moduli"_n] *= (1 - dam) * (1 - dam) + eta;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void
MaterialPhaseField<dim>::computeEffectiveDamageOnQuad(Args && args) {
  // using Mat = Matrix<Real, dim, dim>;

  auto strain = Material::gradUToEpsilon<dim>(args["grad_u"_n]);

  Real trace = strain.trace();
  Real trace_plus = std::max(Real(0.), trace);
  Real trace_minus = std::min(Real(0.), trace);

  Matrix<Real> strain_dev(dim, dim);
  strain_dev = strain - trace / Real(dim) * Matrix<Real>::Identity(dim, dim);

  Real kappa = this->lambda + 2. / Real(dim) * this->mu;

  Real strain_energy_plus = 0.5 * kappa * trace_plus * trace_plus +
                            this->mu * strain_dev.doubleDot(strain_dev);
  Real strain_energy_minus = 0.5 * kappa * trace_minus * trace_minus;

  args["effective_damage"_n] =
      args["damage"_n] * (strain_energy_minus < strain_energy_plus);
}

} // namespace akantu
#endif
