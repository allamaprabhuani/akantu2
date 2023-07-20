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
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialPlastic<dim>::MaterialPlastic(SolidMechanicsModel & model,
                                      const ID & id, const ID & fe_engine_id)
    : MaterialElastic<dim>(model, id, fe_engine_id),
      iso_hardening(this->registerInternal("iso_hardening", 1)),
      inelastic_strain(this->registerInternal("inelastic_strain", dim * dim)),
      plastic_energy(this->registerInternal("plastic_energy", 1)),
      d_plastic_energy(this->registerInternal("d_plastic_energy", 1)) {
  this->registerParam("h", h, Real(0.), _pat_parsable | _pat_modifiable,
                      "Hardening  modulus");
  this->registerParam("sigma_y", sigma_y, Real(0.),
                      _pat_parsable | _pat_modifiable, "Yield stress");

  this->iso_hardening.initializeHistory();
  this->inelastic_strain.initializeHistory();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialPlastic<dim>::getEnergy(const std::string & type) {
  if (type == "plastic") {
    return getPlasticEnergy();
  }
  return MaterialElastic<dim>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <Int dim> Real MaterialPlastic<dim>::getPlasticEnergy() {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();
  Real penergy = 0.;

  for (const auto & type : this->element_filter.elementTypes(dim, _not_ghost)) {
    penergy += fem.integrate(plastic_energy(type, _not_ghost), type, _not_ghost,
                             this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return penergy;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPlastic<dim>::computePotentialEnergy(ElementType el_type) {
  for (auto && [args, epot] :
       zip(getArguments(el_type), this->potential_energy(el_type))) {
    Matrix<Real, dim, dim> elastic_strain =
        args["grad_u"_n] - args["inelastic_strain"_n];

    MaterialElastic<dim>::computePotentialEnergyOnQuad(
        tuple::replace(args, "grad_u"_n = elastic_strain), epot);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialPlastic<dim>::updateEnergies(ElementType el_type) {
  MaterialElastic<dim>::updateEnergies(el_type);

  for (auto && args : zip_append(getArguments(el_type),
                                 "pe"_n = this->plastic_energy(el_type),
                                 "wp"_n = this->d_plastic_energy(el_type))) {

    Matrix<Real, dim, dim> delta_strain =
        args["inelastic_strain"_n] - args["previous_inelastic_strain"_n];

    Matrix<Real, dim, dim> sigma_h = args["sigma"_n] + args["previous_sigma"_n];

    auto && wp = args["wp"_n];
    wp = sigma_h.doubleDot(delta_strain) / 2.;

    args["pe"_n] += wp;
  }
}

/* -------------------------------------------------------------------------- */
template class MaterialPlastic<1>;
template class MaterialPlastic<2>;
template class MaterialPlastic<3>;

} // namespace akantu
