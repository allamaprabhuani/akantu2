/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "local_material_damage.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
LocalMaterialDamage::LocalMaterialDamage(SolidMechanicsModel & model,
                                         const ID & id)
    : Material(model, id), damage(this->registerInternal("damage", 1)) {
  AKANTU_DEBUG_IN();

  this->registerParam("E", E, 0., _pat_parsable, "Young's modulus");
  this->registerParam("nu", nu, 0.5, _pat_parsable, "Poisson's ratio");
  this->registerParam("lambda", lambda, _pat_readable,
                      "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");
  this->registerParam("Yd", Yd, 50., _pat_parsmod);
  this->registerParam("Sd", Sd, 5000., _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LocalMaterialDamage::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  mu = E / (2 * (1 + nu));
  kpa = lambda + 2. / 3. * mu;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LocalMaterialDamage::computeStress(ElementType el_type,
                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  auto dim = this->spatial_dimension;

  for (auto && [grad_u, sigma, dam] :
       zip(make_view(this->getGradU(el_type, ghost_type), dim, dim),
           make_view(this->getStress(el_type, ghost_type), dim, dim),
           damage(el_type, ghost_type))) {
    computeStressOnQuad(grad_u, sigma, dam);
    ++dam;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LocalMaterialDamage::computePotentialEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();
  auto dim = this->spatial_dimension;
  Material::computePotentialEnergy(el_type);

  for (auto && [grad_u, sigma, epot] :
       zip(make_view(this->getGradU(el_type), dim, dim),
           make_view(this->getStress(el_type), dim, dim),
           this->potential_energy(el_type))) {
    computePotentialEnergyOnQuad(grad_u, sigma, epot);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
