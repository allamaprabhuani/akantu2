/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "material_linear_isotropic_hardening.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialLinearIsotropicHardening<dim>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // NOLINTNEXTLINE(bugprone-parent-virtual-call)
  Parent::computeStress(el_type, ghost_type);

  if (this->finite_deformation) { // Finite deformation
    for (auto && [args, green_strain] :
         zip(Parent::getArguments(el_type, ghost_type),
             make_view<dim, dim>((*this->green_strain)(el_type, ghost_type)))) {
      auto & grad_u = args["grad_u"_n];
      auto & previous_grad_u = args["previous_grad_u"_n];

      Material::gradUToE<dim>(grad_u, green_strain);
      auto previous_green_strain = Material::gradUToE<dim>(previous_grad_u);
      auto F_tensor = Material::gradUToF<dim>(grad_u);

      computeStressOnQuad(tuple::append(
          tuple::replace(tuple::replace(args, "grad_u"_n = green_strain),
                         "previous_grad_u"_n = previous_green_strain),
          "F"_n = F_tensor));

      computeStressOnQuad(args);
    }
  } else { // Infinitesimal deformations
    for (auto && args : Parent::getArguments(el_type, ghost_type)) {
      computeStressOnQuad(args);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialLinearIsotropicHardening<spatial_dimension>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && args :
       Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type)) {
    computeTangentModuliOnQuad(args);
  }

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template class MaterialLinearIsotropicHardening<1>;
template class MaterialLinearIsotropicHardening<2>;
template class MaterialLinearIsotropicHardening<3>;
const bool material_is_allocated_plastic_linear_isotropic_hardening
    [[maybe_unused]] = instantiateMaterial<MaterialLinearIsotropicHardening>(
        "plastic_linear_isotropic_hardening");

} // namespace akantu
