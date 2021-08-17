/**
 * @file   material_thermal.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Specialization of the material class for the thermal material
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
#include "material_thermal.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialThermal<dim>::MaterialThermal(SolidMechanicsModel & model,
                                      const ID & id)
    : Material(model, id), delta_T("delta_T", *this),
      sigma_th("sigma_th", *this), use_previous_stress_thermal(false) {
  this->initialize();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialThermal<dim>::MaterialThermal(SolidMechanicsModel & model,
                                      UInt spatial_dimension, const Mesh & mesh,
                                      FEEngine & fe_engine, const ID & id)
    : Material(model, spatial_dimension, mesh, fe_engine, id),
      delta_T("delta_T", *this, dim, fe_engine, this->element_filter),
      sigma_th("sigma_th", *this, dim, fe_engine, this->element_filter),
      use_previous_stress_thermal(false) {
  this->initialize();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialThermal<dim>::initialize() {
  this->registerParam("E", E, Real(0.), _pat_parsable | _pat_modifiable,
                      "Young's modulus");
  this->registerParam("nu", nu, Real(0.5), _pat_parsable | _pat_modifiable,
                      "Poisson's ratio");
  this->registerParam("alpha", alpha, Real(0.), _pat_parsable | _pat_modifiable,
                      "Thermal expansion coefficient");
  this->registerParam("delta_T", delta_T, _pat_parsable | _pat_modifiable,
                      "Uniform temperature field");

  delta_T.initialize(1);
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialThermal<dim>::initMaterial() {
  sigma_th.initialize(1);

  if (use_previous_stress_thermal) {
    sigma_th.initializeHistory();
  }

  Material::initMaterial();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialThermal<dim>::computeStress(ElementType el_type,
                                         GhostType ghost_type) {
  for (auto && tuple : zip(this->delta_T(el_type, ghost_type),
                           this->sigma_th(el_type, ghost_type))) {
    computeStressOnQuad(std::get<1>(tuple), std::get<0>(tuple));
  }
}

/* -------------------------------------------------------------------------- */
INSTANTIATE_MATERIAL_ONLY(MaterialThermal);

} // namespace akantu
