/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_neohookean.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
MaterialNeohookean<spatial_dimension>::MaterialNeohookean(
    SolidMechanicsModel & model, const ID & id)
    : Parent(model, id) {
  this->registerParam("E", E, Real(0.), _pat_parsable | _pat_modifiable,
                      "Young's modulus");
  this->registerParam("nu", nu, Real(0.5), _pat_parsable | _pat_modifiable,
                      "Poisson's ratio");
  this->registerParam("lambda", lambda, _pat_readable,
                      "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");

  this->finite_deformation = true;
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
  PlaneStressToolbox<dim>::computeCauchyStressPlaneStress(el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeCauchyStressPlaneStress(
    ElementType el_type, GhostType ghost_type) {
  for (auto && [grad_u, piola, sigma, c33] :
       zip(make_view<2, 2>((*this->gradu)(el_type, ghost_type)),
           make_view<2, 2>((*this->piola_kirchhoff_2)(el_type, ghost_type)),
           make_view<2, 2>((*this->stress)(el_type, ghost_type)),
           make_view((*this->third_axis_deformation)(el_type, ghost_type)))) {
    StoCauchy<2>(gradUToF<2>(grad_u), piola, sigma, c33);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialNeohookean<dim>::computeStress(ElementType el_type,
                                            GhostType ghost_type) {
  PlaneStressToolbox<dim>::computeStress(el_type, ghost_type);
  for (auto && args : this->getArguments(el_type, ghost_type)) {
    computeStressOnQuad(args);
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

  for (auto && args : getArguments(el_type, ghost_type)) {
    computeThirdAxisDeformationOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialNeohookean<dim>::computePotentialEnergy(ElementType el_type) {
  Parent::computePotentialEnergy(el_type);
  auto && arguments = Parent::getArguments(el_type, _not_ghost);

  for (auto && [args, epot] :
       zip(arguments, (*this->potential_energy)(el_type, _not_ghost))) {
    this->computePotentialEnergyOnQuad(args, epot);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialNeohookean<dim>::computeTangentModuli(ElementType el_type,
                                                   Array<Real> & tangent_matrix,
                                                   GhostType ghost_type) {
  auto && arguments =
      Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type);
  Parent::computeStress(el_type, ghost_type);

  for (auto && args : arguments) {
    computeTangentModuliOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialNeohookean<dim>::getPushWaveSpeed(
    const Element & /*element*/) const {
  return std::sqrt((this->lambda + 2 * this->mu) / this->rho);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialNeohookean<dim>::getShearWaveSpeed(
    const Element & /*element*/) const {
  return std::sqrt(this->mu / this->rho);
}

/* -------------------------------------------------------------------------- */
template class MaterialNeohookean<1>;
template class MaterialNeohookean<2>;
template class MaterialNeohookean<3>;
static bool material_is_allocated_neohookean =
    instantiateMaterial<MaterialNeohookean>("neohookean");

} // namespace akantu
