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
#include "material_elastic.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialElastic<dim>::MaterialElastic(SolidMechanicsModel & model,
                                      const ID & id, const ID & fe_engine_id)
    : Parent(model, id, fe_engine_id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialElastic<dim>::initialize() {
  this->registerParam("lambda", lambda, _pat_readable,
                      "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialElastic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  Parent::initMaterial();

  if (dim == 1) {
    this->nu = 0.;
  }

  this->updateInternalParameters();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialElastic<dim>::updateInternalParameters() {
  MaterialThermal<dim>::updateInternalParameters();

  this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  this->mu = this->E / (2 * (1 + this->nu));

  this->kpa = this->lambda + 2. / 3. * this->mu;

  this->was_stiffness_assembled = false;
}

/* -------------------------------------------------------------------------- */
template <> void MaterialElastic<2>::updateInternalParameters() {
  MaterialThermal<2>::updateInternalParameters();

  this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  this->mu = this->E / (2 * (1 + this->nu));

  if (this->plane_stress) {
    this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - this->nu));
  }

  this->kpa = this->lambda + 2. / 3. * this->mu;

  this->was_stiffness_assembled = false;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElastic<dim>::computeStress(ElementType el_type,
                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent::computeStress(el_type, ghost_type);

  auto && arguments = Parent::getArguments(el_type, ghost_type);

  if (not this->finite_deformation) {
    for (auto && args : arguments) {
      this->computeStressOnQuad(args);
    }
  } else {
    for (auto && args : arguments) {
      auto && E = this->template gradUToE<dim>(args["grad_u"_n]);
      this->computeStressOnQuad(tuple::replace(args, "grad_u"_n = E));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElastic<dim>::computeTangentModuli(ElementType el_type,
                                                Array<Real> & tangent_matrix,
                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto && arguments =
      Parent::getArgumentsTangent(tangent_matrix, el_type, ghost_type);

  for (auto && args : arguments) {
    this->computeTangentModuliOnQuad(args);
  }

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialElastic<dim>::getPushWaveSpeed(const Element &) const {
  return sqrt((lambda + 2 * mu) / this->rho);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialElastic<dim>::getShearWaveSpeed(const Element &) const {
  return sqrt(mu / this->rho);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElastic<dim>::computePotentialEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();

  // needs to be implemented
  // MaterialThermal<dim>::computePotentialEnergy(el_type);
  auto && arguments = Parent::getArguments(el_type, _not_ghost);

  if (not this->finite_deformation) {
    for (auto && [args, epot] :
         zip(arguments, (*this->potential_energy)(el_type, _not_ghost))) {
      this->computePotentialEnergyOnQuad(args, epot);
    }
  } else {
    for (auto && [args, epot] :
         zip(arguments, (*this->potential_energy)(el_type, _not_ghost))) {
      auto && E = this->template gradUToE<dim>(args["grad_u"_n]);
      this->computePotentialEnergyOnQuad(tuple::replace(args, "grad_u"_n = E),
                                         epot);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElastic<dim>::computePotentialEnergyByElement(
    const Element & element, Vector<Real> & epot_on_quad_points) {
  auto type = element.type;
  auto gradu_view = make_view<dim, dim>((*this->gradu)(type));
  auto stress_view = make_view<dim, dim>((*this->stress)(type));

  if (this->finite_deformation) {
    stress_view = make_view<dim, dim>((*this->piola_kirchhoff_2)(type));
  }

  auto nb_quadrature_points = this->getFEEngine().getNbIntegrationPoints(type);

  auto gradu_it = gradu_view.begin() + element.element * nb_quadrature_points;
  auto gradu_end = gradu_it + nb_quadrature_points;
  auto stress_it = stress_view.begin() + element.element * nb_quadrature_points;
  auto stress_end = stress_it + nb_quadrature_points;

  auto epot_quad = epot_on_quad_points.begin();

  if (this->finite_deformation) {
    for (auto && data : zip("grad_u"_n = range(gradu_it, gradu_end),
                        "sigma"_n = range(stress_it, stress_end),
                        "Epot"_n = epot_on_quad_points)) {
      auto E = this->template gradUToE<dim>(data["grad_u"_n]);
      this->computePotentialEnergyOnQuad(tuple::replace(data, "grad_u"_n = E),
                                         data["Epot"_n]);
    }
  } else {
    for (auto && data : zip("grad_u"_n = range(gradu_it, gradu_end),
                        "sigma"_n = range(stress_it, stress_end),
                        "Epot"_n = epot_on_quad_points)) {
      this->computePotentialEnergyOnQuad(data, data["Epot"_n]);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
Real MaterialElastic<1>::getPushWaveSpeed(const Element & /*element*/) const {
  return std::sqrt(this->E / this->rho);
}

template <>
Real MaterialElastic<1>::getShearWaveSpeed(const Element & /*element*/) const {
  AKANTU_EXCEPTION("There is no shear wave speed in 1D");
}

/* -------------------------------------------------------------------------- */
template class MaterialElastic<1>;
template class MaterialElastic<2>;
template class MaterialElastic<3>;

static bool material_is_allocated_elastic =
    instantiateMaterial<MaterialElastic>("elastic");

} // namespace akantu
