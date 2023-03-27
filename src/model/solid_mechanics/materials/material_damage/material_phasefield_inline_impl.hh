/**
 * @file   material_phasefield_inline_impl.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Fri Apr 02 2021
 *
 * @brief  Implementation of the inline functions of the material phasefield
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

#include "material_phasefield.hh"
#include <algorithm>

#ifndef AKANTU_MATERIAL_PHASEFIELD_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_PHASEFIELD_INLINE_IMPL_HH_
/* -------------------------------------------------------------------------- */
namespace akantu {
template <Int dim>
template <class Args>
inline void MaterialPhaseField<dim>::computeStressOnQuad(Args && args) {
  MaterialElastic<dim>::computeStressOnQuad(args);

  auto && dam = args["effective_damage"_n];
  args["sigma"_n] *= (1 - dam) * (1 - dam) + eta;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
void MaterialPhaseField<dim>::computeTangentModuliOnQuad(Args && args) {
  MaterialElastic<dim>::computeTangentModuliOnQuad(args);

  auto dam = args["effective_damage"_n];
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
  strain_dev = strain - trace / static_cast<Real>(dim) *
                            Matrix<Real>::Identity(dim, dim);

  Real Kn = this->lambda + 2. * this->mu / static_cast<Real>(dim);

  Real strain_energy_plus = 0.5 * Kn * trace_plus * trace_plus +
                            this->mu * strain_dev.doubleDot(strain_dev);
  Real strain_energy_minus = 0.5 * Kn * trace_minus * trace_minus;

  // Mat strain_dir;
  // Vector<Real, dim> strain_values;
  // strain.eig(strain_values, strain_dir);

  // Mat strain_diag_plus;
  // Mat strain_diag_minus;

  // strain_diag_plus.zero();
  // strain_diag_minus.zero();

  // for (UInt i = 0; i < dim; i++) {
  //   strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
  //   strain_diag_minus(i, i) = std::min(Real(0.), strain_values(i));
  // }

  // Mat strain_plus = strain_dir * strain_diag_plus * strain_dir.transpose();
  // Mat strain_minus = strain_dir * strain_diag_minus * strain_dir.transpose();

  // auto trace_plus = std::max(Real(0.), strain.trace());
  // auto trace_minus = std::min(Real(0.), strain.trace());

  // Mat sigma_plus =
  //     Mat::Identity() * trace_plus * this->lambda + 2. * this->mu *
  //     strain_plus;
  // Mat sigma_minus = Mat::Identity() * trace_minus * this->lambda +
  //                   2. * this->mu * strain_minus;

  // auto strain_energy_plus = sigma_plus.doubleDot(strain_plus) / 2.;
  // auto strain_energy_minus = sigma_minus.doubleDot(strain_minus) / 2.;

  args["effective_damage"_n] =
      args["damage"_n] * (strain_energy_minus < strain_energy_plus);
}

} // namespace akantu
#endif
