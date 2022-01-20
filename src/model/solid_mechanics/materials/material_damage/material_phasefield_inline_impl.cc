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

/* -------------------------------------------------------------------------- */
namespace akantu {
template <Int dim>
template <class Args>
inline void MaterialPhaseField<dim>::computeStressOnQuad(Args && args) {
  MaterialElastic<dim>::computeStressOnQuad(args);

  auto && dam = tuple::get<"damage"_h>(args);
  tuple::get<"sigma"_h>(args) *= (1 - dam) * (1 - dam) + eta;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
void MaterialPhaseField<dim>::computeTangentModuliOnQuad(Args && args) {
  MaterialElastic<dim>::computeTangentModuliOnQuad(args);

  auto dam = tuple::get<"damage"_h>(args);
  tuple::get<"tangent_moduli"_h>(args) *= (1 - dam) * (1 - dam) + eta;
}

} // namespace akantu
