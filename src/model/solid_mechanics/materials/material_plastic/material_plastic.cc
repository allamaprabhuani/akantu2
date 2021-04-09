/**
 * @file   material_plastic.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Apr 07 2014
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Implemantation of the akantu::MaterialPlastic class
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
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialPlastic<spatial_dimension>::MaterialPlastic(SolidMechanicsModel & model,
                                                    const ID & id,
                                                    const ID & fe_engine_id)
    : MaterialElastic<spatial_dimension>(model, id, fe_engine_id),
      iso_hardening("iso_hardening", *this, fe_engine_id, this->element_filter),
      inelastic_strain("inelastic_strain", *this, fe_engine_id,
                       this->element_filter),
      plastic_energy("plastic_energy", *this, fe_engine_id,
                     this->element_filter),
      d_plastic_energy("d_plastic_energy", *this, fe_engine_id,
                       this->element_filter) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::initialize() {
  this->registerParam("h", h, Real(0.), _pat_parsable | _pat_modifiable,
                      "Hardening  modulus");
  this->registerParam("sigma_y", sigma_y, Real(0.),
                      _pat_parsable | _pat_modifiable, "Yield stress");

  this->iso_hardening.initialize(1);
  this->iso_hardening.initializeHistory();

  this->plastic_energy.initialize(1);
  this->d_plastic_energy.initialize(1);

  this->use_previous_stress = true;
  this->use_previous_gradu = true;
  this->use_previous_stress_thermal = true;

  this->inelastic_strain.initialize(spatial_dimension * spatial_dimension);
  this->inelastic_strain.initializeHistory();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getEnergy(const std::string & type) {
  if (type == "plastic") {
    return getPlasticEnergy();
  }
  return MaterialElastic<spatial_dimension>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getPlasticEnergy() {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();
  Real penergy = 0.;

  for (auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    penergy += fem.integrate(plastic_energy(type, _not_ghost), type, _not_ghost,
                             this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return penergy;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  auto epot = this->potential_energy(el_type).begin();

  auto inelastic_strain_it = this->inelastic_strain(el_type).begin(
      spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  Matrix<Real> elastic_strain(spatial_dimension, spatial_dimension);
  elastic_strain.copy(grad_u);
  elastic_strain -= *inelastic_strain_it;

  MaterialElastic<spatial_dimension>::computePotentialEnergyOnQuad(
      elastic_strain, sigma, *epot);

  ++epot;
  ++inelastic_strain_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::updateEnergies(ElementType el_type) {
  AKANTU_DEBUG_IN();

  MaterialElastic<spatial_dimension>::updateEnergies(el_type);

  auto pe_it = this->plastic_energy(el_type).begin();
  auto wp_it = this->d_plastic_energy(el_type).begin();

  auto inelastic_strain_it = this->inelastic_strain(el_type).begin(
      spatial_dimension, spatial_dimension);

  auto previous_inelastic_strain_it =
      this->inelastic_strain.previous(el_type).begin(spatial_dimension,
                                                     spatial_dimension);

  auto previous_sigma = this->stress.previous(el_type).begin(spatial_dimension,
                                                             spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  Matrix<Real> delta_strain_it(*inelastic_strain_it);
  delta_strain_it -= *previous_inelastic_strain_it;

  Matrix<Real> sigma_h(sigma);
  sigma_h += *previous_sigma;

  *wp_it = .5 * sigma_h.doubleDot(delta_strain_it);

  *pe_it += *wp_it;

  ++pe_it;
  ++wp_it;
  ++inelastic_strain_it;
  ++previous_inelastic_strain_it;
  ++previous_sigma;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL_ONLY(MaterialPlastic);

} // namespace akantu
