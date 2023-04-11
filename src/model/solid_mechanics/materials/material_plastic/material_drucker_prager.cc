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
#include "material_drucker_prager.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template <Int dim>
MaterialDruckerPrager<spatial_dimension>::MaterialDruckerPrager(
    SolidMechanicsModel & model, const ID & id, const ID & fe_engine_id)
    : MaterialPlastic<spatial_dimension>(model, id, fe_engine_id) {
  this->initialize();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialDruckerPrager<dim>::initialize() {
  this->registerParam("phi", phi, Real(0.), _pat_parsable | _pat_modifiable,
                      "Internal friction angle in degrees");
  this->registerParam("fc", fc, Real(1.), _pat_parsable | _pat_modifiable,
                      "Compressive strength");
  this->registerParam("radial_return", radial_return_mapping, bool(true),
                      _pat_parsable | _pat_modifiable, "Radial return mapping");
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialDruckerPrager<dim>::updateInternalParameters() {
  Parent::updateInternalParameters();

  // compute alpha and k parameters for Drucker-Prager
  Real phi_radian = this->phi * M_PI / 180.;
  this->alpha = (6. * sin(phi_radian)) / (3. - sin(phi_radian));
  Real cohesion = this->fc * (1. - sin(phi_radian)) / (2. * cos(phi_radian));
  this->k = (6. * cohesion * cos(phi_radian)) / (3. - sin(phi_radian));
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialDruckerPrager<dim>::computeStress(ElementType el_type,
                                               GhostType ghost_type) {
  if (this->finite_deformation) {
    AKANTU_TO_IMPLEMENT();
  }

  MaterialThermal<dim>::computeStress(el_type, ghost_type);

  // infinitesimal and finite deformation
  for (auto && args : Parent::getArguments(el_type, ghost_type)) {
    computeStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template class MaterialDruckerPrager<1>;
template class MaterialDruckerPrager<2>;
template class MaterialDruckerPrager<3>;
static bool material_is_allocated_plastic_drucker_prager =
    instantiateMaterial<MaterialDruckerPrager>("plastic_drucker_prager");

} // namespace akantu
