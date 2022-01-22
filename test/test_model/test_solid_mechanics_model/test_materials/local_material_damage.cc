/**
 * @file   local_material_damage.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Oct 19 2014
 * @date last modification:  Fri May 03 2019
 *
 * @brief  Specialization of the material class for the damage material
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
#include "local_material_damage.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
LocalMaterialDamage::LocalMaterialDamage(SolidMechanicsModel & model,
                                         const ID & id)
    : Material(model, id), damage("damage", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("E", E, 0., _pat_parsable, "Young's modulus");
  this->registerParam("nu", nu, 0.5, _pat_parsable, "Poisson's ratio");
  this->registerParam("lambda", lambda, _pat_readable,
                      "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");
  this->registerParam("Yd", Yd, 50., _pat_parsmod);
  this->registerParam("Sd", Sd, 5000., _pat_parsmod);

  damage.initialize(1);

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

  for (auto && data :
       zip(make_view(this->gradu(el_type, ghost_type), dim, dim),
           make_view(this->stress(el_type, ghost_type), dim, dim),
           damage(el_type, ghost_type))) {
    auto && grad_u = std::get<0>(data);
    auto && sigma = std::get<1>(data);
    auto && dam = std::get<2>(data);

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

  for (auto && data : zip(make_view(this->gradu(el_type), dim, dim),
                          make_view(this->stress(el_type), dim, dim),
                          potential_energy(el_type))) {
    auto && grad_u = std::get<0>(data);
    auto && sigma = std::get<1>(data);
    auto && epot = std::get<2>(data);
    computePotentialEnergyOnQuad(grad_u, sigma, epot);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
