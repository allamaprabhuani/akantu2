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
template <Int spatial_dimension>
MaterialPlastic<spatial_dimension>::MaterialPlastic(SolidMechanicsModel & model,
                                                    const ID & id)
    : MaterialElastic<spatial_dimension>(model, id),
      iso_hardening("iso_hardening", *this),
      inelastic_strain("inelastic_strain", *this),
      plastic_energy("plastic_energy", *this),
      d_plastic_energy("d_plastic_energy", *this) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

template <Int spatial_dimension>
MaterialPlastic<spatial_dimension>::MaterialPlastic(SolidMechanicsModel & model,
                                                    UInt dim, const Mesh & mesh,
                                                    FEEngine & fe_engine,
                                                    const ID & id)
    : MaterialElastic<spatial_dimension>(model, dim, mesh, fe_engine, id),
      iso_hardening("iso_hardening", *this, dim, fe_engine,
                    this->element_filter),
      inelastic_strain("inelastic_strain", *this, dim, fe_engine,
                       this->element_filter),
      plastic_energy("plastic_energy", *this, dim, fe_engine,
                     this->element_filter),
      d_plastic_energy("d_plastic_energy", *this, dim, fe_engine,
                       this->element_filter) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
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

  this->inelastic_strain.initialize(spatial_dimension * spatial_dimension);
  this->inelastic_strain.initializeHistory();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getEnergy(const std::string & type) {
  if (type == "plastic") {
    return getPlasticEnergy();
  }
  return MaterialElastic<spatial_dimension>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getPlasticEnergy() {
  AKANTU_DEBUG_IN();

  Real penergy = 0.;

  for (const auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    penergy +=
        this->fem.integrate(plastic_energy(type, _not_ghost), type, _not_ghost,
                            this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return penergy;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPlastic<dim>::computePotentialEnergy(ElementType el_type) {
  auto epot = this->potential_energy(el_type).begin();

  for (auto && args : getArguments(el_type)) {
    Matrix<Real, dim, dim> elastic_strain =
        args["grad_u"_n] - args["inelastic_strain"_n];

    MaterialElastic<dim>::computePotentialEnergyOnQuad(
        tuple::replace(args, "grad_u"_n = elastic_strain), *epot);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPlastic<dim>::updateEnergies(ElementType el_type) {
  MaterialElastic<dim>::updateEnergies(el_type);

  for (auto && args : zip_append(getArguments(el_type),
                                 "pe"_n = this->plastic_energy(el_type),
                                 "wp"_n = this->d_plastic_energy(el_type))) {

    Matrix<Real, dim, dim> delta_strain_it =
        args["inelastic_strain"_n] - args["previous_inelastic_strain"_n];

    Matrix<Real, dim, dim> sigma_h = args["sigma"_n] + args["previous_sigma"_n];

    auto && wp = args["wp"_n];
    wp = .5 * sigma_h.doubleDot(delta_strain_it);

    args["pe"_n] += wp;
  }
}

/* -------------------------------------------------------------------------- */
template class MaterialPlastic<1>;
template class MaterialPlastic<2>;
template class MaterialPlastic<3>;

} // namespace akantu
