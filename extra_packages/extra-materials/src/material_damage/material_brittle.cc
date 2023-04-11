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
#include "material_brittle.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
MaterialBrittle<spatial_dimension>::MaterialBrittle(SolidMechanicsModel & model,
                                                    const ID & id)
    : MaterialDamage<spatial_dimension>(model, id),
      strain_rate_brittle("strain_rate_brittle", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("S_0", S_0, 157e6, _pat_parsable | _pat_modifiable);
  this->registerParam("E_0", E_0, 27e3, _pat_parsable, "Strain rate threshold");
  this->registerParam("A", A, 1.622e-5, _pat_parsable,
                      "Polynome cubic constant");
  this->registerParam("B", B, -1.3274, _pat_parsable,
                      "Polynome quadratic constant");
  this->registerParam("C", C, 3.6544e4, _pat_parsable,
                      "Polynome linear constant");
  this->registerParam("D", D, -181.38e6, _pat_parsable, "Polynome constant");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialBrittle<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();

  this->strain_rate_brittle.initialize(spatial_dimension * spatial_dimension);

  updateInternalParameters();

  // this->Yd.resize();
  // const Mesh & mesh = this->model.getFEEngine().getMesh();

  // Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  // Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);

  // for(; it != last_type; ++it) {
  //   UInt nb_element  = this->element_filter(*it).getSize();
  //   UInt nb_quad = this->model.getFEEngine().getNbQuadraturePoints(*it);

  //   Array <Real> & Yd_rand_vec = Yd_rand(*it);
  //   for(UInt e = 0; e < nb_element; ++e) {
  //     Real rand_part = (2 * drand48()-1) * Yd_randomness * Yd;

  //     for(UInt q = 0; q < nb_quad; ++q)
  //    Yd_rand_vec(nb_quad*e+q,0) = Yd + rand_part;
  //   }
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialBrittle<spatial_dimension>::updateInternalParameters() {
  MaterialDamage<spatial_dimension>::updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialBrittle<spatial_dimension>::computeStress(ElementType el_type,
                                                       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).data();

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

  Real sigma_equivalent = 0;
  Real fracture_stress = 0;
  Matrix<Real> & grad_v = *strain_rate_it;

  computeStressOnQuad(grad_u, grad_v, sigma, *dam, sigma_equivalent,
                      fracture_stress);

  ++strain_rate_it;
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(brittle, MaterialBrittle);

} // namespace akantu
