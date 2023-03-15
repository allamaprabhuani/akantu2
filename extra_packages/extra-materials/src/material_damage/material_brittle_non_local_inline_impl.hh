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
MaterialBrittleNonLocal<spatial_dimension>::MaterialBrittleNonLocal(
    SolidMechanicsModel & model, const ID & id)
    : Material(model, id), MaterialBrittleNonLocalParent(model, id),
      Sigma_max("Sigma max", *this), Sigma_maxnl("Sigma_max non local", *this),
      Sigma_fracture("Sigma_fracture", *this) {
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->Sigma_max.initialize(1);
  this->Sigma_maxnl.initialize(1);
  this->Sigma_fracture.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialBrittleNonLocal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  this->model.getNonLocalManager().registerNonLocalVariable(
      this->Sigma_max.getName(), Sigma_maxnl.getName(), 1);
  MaterialBrittleNonLocalParent::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialBrittleNonLocal<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).data();
  Real * Sigma_maxt = this->Sigma_max(el_type, ghost_type).data();
  Real * fracture_stress = this->Sigma_fracture(el_type, ghost_type).data();

  Array<Real> & velocity = this->model.getVelocity();
  Array<Real> & strain_rate_brittle =
      this->strain_rate_brittle(el_type, ghost_type);
  Array<Idx> & elem_filter = this->element_filter(el_type, ghost_type);

  this->model.getFEEngine().gradientOnIntegrationPoints(
      velocity, strain_rate_brittle, spatial_dimension, el_type, ghost_type,
      elem_filter);

  Array<Real>::iterator<Matrix<Real>> strain_rate_it =
      this->strain_rate_brittle(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & grad_v = *strain_rate_it;

  MaterialBrittle<spatial_dimension>::computeStressOnQuad(
      grad_u, grad_v, sigma, *dam, *Sigma_maxt, *fracture_stress);
  ++dam;
  ++strain_rate_it;
  ++Sigma_maxt;
  ++fracture_stress;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialBrittleNonLocal<spatial_dimension>::computeNonLocalStress(
    ElementType type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(type, ghost_type).data();
  Real * Sigma_maxnlt = this->Sigma_maxnl(type, ghost_type).data();
  Real * fracture_stress = this->Sigma_fracture(type, ghost_type).data();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(type, ghost_type);
  this->computeDamageAndStressOnQuad(sigma, *dam, *Sigma_maxnlt,
                                     *fracture_stress);

  ++dam;
  ++Sigma_maxnlt;
  ++fracture_stress;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialBrittleNonLocal<
    spatial_dimension>::nonLocalVariableToNeighborhood() {
  this->model.getNonLocalManager().nonLocalVariableToNeighborhood(
      Sigma_maxnl.getName(), this->name);
}
