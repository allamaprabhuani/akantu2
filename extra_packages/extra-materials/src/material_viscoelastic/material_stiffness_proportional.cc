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
#include "material_stiffness_proportional.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
MaterialStiffnessProportional<spatial_dimension>::MaterialStiffnessProportional(
    SolidMechanicsModel & model, const ID & id)
    : MaterialElastic<spatial_dimension>(model, id),
      stress_viscosity("stress_viscosity", *this),
      stress_elastic("stress_elastic", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("Alpha", alpha, 0., _pat_parsable | _pat_modifiable,
                      "Artificial viscous ratio");

  this->stress_viscosity.initialize(spatial_dimension * spatial_dimension);
  this->stress_elastic.initialize(spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialStiffnessProportional<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic<spatial_dimension>::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialStiffnessProportional<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  Array<Idx> & elem_filter = this->element_filter(el_type, ghost_type);
  Array<Real> & stress_visc = stress_viscosity(el_type, ghost_type);
  Array<Real> & stress_el = stress_elastic(el_type, ghost_type);

  MaterialElastic<spatial_dimension>::computeStress(el_type, ghost_type);

  Array<Real> & velocity = this->model.getVelocity();

  Array<Real> strain_rate(0, spatial_dimension * spatial_dimension,
                          "strain_rate");

  this->model.getFEEngine().gradientOnIntegrationPoints(
      velocity, strain_rate, spatial_dimension, el_type, ghost_type,
      elem_filter);

  Array<Real>::matrix_iterator strain_rate_it =
      strain_rate.begin(spatial_dimension, spatial_dimension);
  Array<Real>::matrix_iterator stress_visc_it =
      stress_visc.begin(spatial_dimension, spatial_dimension);
  Array<Real>::matrix_iterator stress_el_it =
      stress_el.begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  Matrix<Real> & grad_v = *strain_rate_it;
  Matrix<Real> & sigma_visc = *stress_visc_it;
  Matrix<Real> & sigma_el = *stress_el_it;

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_v, sigma_visc);

  sigma_visc *= alpha;
  sigma_el.copy(sigma);
  sigma += sigma_visc;

  ++strain_rate_it;
  ++stress_visc_it;
  ++stress_el_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialStiffnessProportional<spatial_dimension>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  Array<Real> & stress_el = stress_elastic(el_type);
  Array<Real>::matrix_iterator stress_el_it =
      stress_el.begin(spatial_dimension, spatial_dimension);

  Real * epot = this->potential_energy(el_type).data();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);
  Matrix<Real> & sigma_el = *stress_el_it;
  MaterialElastic<spatial_dimension>::computePotentialEnergyOnQuad(
      grad_u, sigma_el, *epot);
  epot++;
  ++stress_el_it;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
INSTANTIATE_MATERIAL(ve_stiffness_prop, MaterialStiffnessProportional);

} // namespace akantu
