/**
 * @file   material_linear_isotropic_hardening.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Benjamin Paccaud <benjamin.paccaud@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Apr 07 2014
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Specialization of the material class for isotropic finite deformation
 * linear hardening plasticity
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_linear_isotropic_hardening.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialLinearIsotropicHardening<dim>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // NOLINTNEXTLINE(bugprone-parent-virtual-call)
  MaterialThermal<dim>::computeStress(el_type, ghost_type);

  if (this->finite_deformation) { // Finite deformation
    for (auto && data :
         zip(Parent::getArguments(el_type, ghost_type),
             make_view<dim, dim>(this->green_strain(el_type, ghost_type)))) {
      auto && args = std::get<0>(data);
      auto & green_strain = std::get<1>(data);
      auto & grad_u = tuple::get<"grad_u"_h>(args);
      auto & previous_grad_u = tuple::get<"previous_grad_u"_h>(args);

      Material::gradUToE<dim>(grad_u, green_strain);
      Matrix<Real, dim, dim> previous_green_strain =
          Material::gradUToE<dim>(previous_grad_u);
      Matrix<Real, dim, dim> F_tensor = Material::gradUToF<dim>(grad_u);

      computeStressOnQuad(
          tuple::append(tuple::replace<"previous_grad_u"_h>(
                            tuple::replace<"grad_u"_h>(args, green_strain),
                            previous_green_strain),
                        tuple::get<"F"_h>() = F_tensor));

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
static bool material_is_allocated_plastic_linear_isotropic_hardening =
    instantiateMaterial<MaterialLinearIsotropicHardening>(
        "plastic_linear_isotropic_hardening");

} // namespace akantu
