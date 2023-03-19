/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
} // namespace akantu

#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#include <string>
#endif

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
MaterialDamageIterativeNonLocal<spatial_dimension>::
    MaterialDamageIterativeNonLocal(SolidMechanicsModel & model, const ID & id)
    : MaterialDamageIterativeNonLocalParent(model, id),
      grad_u_nl("grad_u non local", *this) {
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->grad_u_nl.initialize(spatial_dimension * spatial_dimension);
  this->model.getNonLocalManager().registerNonLocalVariable(
      this->gradu.getName(), grad_u_nl.getName(),
      spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialDamageIterativeNonLocal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamageIterativeNonLocalParent::initMaterial();

  this->model.getNonLocalManager().nonLocalVariableToNeighborhood(
      grad_u_nl.getName(), this->name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialDamageIterativeNonLocal<spatial_dimension>::computeStress(
    ElementType /*type*/, GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialDamageIterativeNonLocal<spatial_dimension>::computeNonLocalStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// compute the stress (based on the elastic law)
  MaterialDamage<spatial_dimension>::computeStress(el_type, ghost_type);

  /// multiply the stress by (1-d) to get the effective stress
  Real * dam = this->damage(el_type, ghost_type).data();
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeDamageAndStressOnQuad(sigma, *dam);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  /// compute the normalized equivalent stress
  this->computeNormalizedEquivalentStress(this->grad_u_nl(el_type, ghost_type),
                                          el_type, ghost_type);
  /// find the maximum
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}
