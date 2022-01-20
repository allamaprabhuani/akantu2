/**
 * @file   material_neohookean.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Apr 08 2013
 * @date last modification: Thu Feb 20 2020
 *
 * @brief  Specialization of the material class for finite deformation
 * neo-hookean material
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
#include "material_neohookean.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
MaterialNeohookean<spatial_dimension>::MaterialNeohookean(
    SolidMechanicsModel & model, const ID & id)
    : Parent(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("E", E, Real(0.), _pat_parsable | _pat_modifiable,
                      "Young's modulus");
  this->registerParam("nu", nu, Real(0.5), _pat_parsable | _pat_modifiable,
                      "Poisson's ratio");
  this->registerParam("lambda", lambda, _pat_readable,
                      "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");

  this->finite_deformation = true;
  this->initialize_third_axis_deformation = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialNeohookean<spatial_dimension>::initMaterial() {
  PlaneStressToolbox<spatial_dimension>::initMaterial();
  if (spatial_dimension == 1) {
    nu = 0.;
  }
  this->updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template <> void MaterialNeohookean<2>::initMaterial() {
  PlaneStressToolbox<2>::initMaterial();
  this->updateInternalParameters();

  if (this->plane_stress) {
    this->third_axis_deformation.setDefaultValue(1.);
  }
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialNeohookean<spatial_dimension>::updateInternalParameters() {
  lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  mu = E / (2 * (1 + nu));

  kpa = lambda + 2. / 3. * mu;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialNeohookean<dim>::computeCauchyStressPlaneStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  PlaneStressToolbox<dim>::computeCauchyStressPlaneStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeCauchyStressPlaneStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto gradu_it = this->gradu(el_type, ghost_type).begin(2, 2);
  auto gradu_end = this->gradu(el_type, ghost_type).end(2, 2);
  auto piola_it = this->piola_kirchhoff_2(el_type, ghost_type).begin(2, 2);
  auto stress_it = this->stress(el_type, ghost_type).begin(2, 2);
  auto c33_it = this->third_axis_deformation(el_type, ghost_type).begin();

  for (; gradu_it != gradu_end; ++gradu_it, ++piola_it, ++stress_it, ++c33_it) {
    auto && grad_u = *gradu_it;
    auto && piola = *piola_it;
    auto && sigma = *stress_it;

    StoCauchy<2>(gradUToF<2>(grad_u), piola, sigma, *c33_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialNeohookean<dim>::computeStress(ElementType el_type,
                                            GhostType ghost_type) {
  for (auto && args : this->getArguments(el_type, ghost_type)) {
    computeStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeStress(ElementType el_type,
                                          GhostType ghost_type) {
  auto && arguments = getArguments(el_type, ghost_type);
  if (this->plane_stress) {
    PlaneStressToolbox<2>::computeStress(el_type, ghost_type);

    for (auto && args : zip_replace<"C33"_h>(
             arguments, this->third_axis_deformation(el_type, ghost_type))) {
      computeStressOnQuad(args);
    }
  } else {
    for (auto && args : arguments) {
      computeStressOnQuad(args);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialNeohookean<dim>::computeThirdAxisDeformation(
    ElementType /*el_type*/, GhostType /*ghost_type*/) {}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeThirdAxisDeformation(ElementType el_type,
                                                        GhostType ghost_type) {
  AKANTU_DEBUG_ASSERT(this->plane_stress, "The third component of the strain "
                                          "can only be computed for 2D "
                                          "problems in Plane Stress!!");

  auto && arguments = getArguments(el_type, ghost_type);
  for (auto && args : zip_replace<"C33"_h>(
           arguments, this->third_axis_deformation(el_type, ghost_type))) {
    computeThirdAxisDeformationOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();
  constexpr Int dim = 2;

  Material::computePotentialEnergy(el_type);

  Array<Real>::scalar_iterator epot = this->potential_energy(el_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  computePotentialEnergyOnQuad(grad_u, *epot);
  ++epot;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  auto && arguments = getArgumentsTangent(tangent_matrix, el_type, ghost_type);

  for (auto && args : arguments) {
    computeTangentModuliOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeTangentModuli(ElementType el_type,
                                                 Array<Real> & tangent_matrix,
                                                 GhostType ghost_type) {
  auto && arguments = getArgumentsTangent(tangent_matrix, el_type, ghost_type);

  if (this->plane_stress) {
    PlaneStressToolbox<2>::computeStress(el_type, ghost_type);

    for (auto && args : zip_replace<"C33"_h>(
             arguments, this->third_axis_deformation(el_type, ghost_type))) {
      computeTangentModuliOnQuad(args);
    }
  } else {
    for (auto && args : arguments) {
      computeTangentModuliOnQuad(args);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
Real MaterialNeohookean<spatial_dimension>::getPushWaveSpeed(
    const Element & /*element*/) const {
  return sqrt((this->lambda + 2 * this->mu) / this->rho);
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
Real MaterialNeohookean<spatial_dimension>::getShearWaveSpeed(
    const Element & /*element*/) const {
  return sqrt(this->mu / this->rho);
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(neohookean, MaterialNeohookean);

} // namespace akantu
