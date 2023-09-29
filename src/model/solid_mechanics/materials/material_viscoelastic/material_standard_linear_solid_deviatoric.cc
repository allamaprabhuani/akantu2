/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_standard_linear_solid_deviatoric.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialStandardLinearSolidDeviatoric<
    dim>::MaterialStandardLinearSolidDeviatoric(SolidMechanicsModel & model,
                                                const ID & id)
    : MaterialElastic<dim>(model, id),
      stress_dev(this->registerInternal("stress_dev", dim * dim)),
      history_integral(this->registerInternal("history_integral", dim * dim)),
      dissipated_energy(this->registerInternal("dissipated_energy", 1)) {
  this->registerParam("Eta", eta, Real(1.), _pat_parsable | _pat_modifiable,
                      "Viscosity");
  this->registerParam("Ev", Ev, Real(1.), _pat_parsable | _pat_modifiable,
                      "Stiffness of the viscous element");
  this->registerParam("Einf", E_inf, Real(1.), _pat_readable,
                      "Stiffness of the elastic element");
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialStandardLinearSolidDeviatoric<dim>::updateInternalParameters() {
  MaterialElastic<dim>::updateInternalParameters();
  E_inf = this->E - this->Ev;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialStandardLinearSolidDeviatoric<dim>::setToSteadyState(
    ElementType el_type, GhostType ghost_type) {
  /// Loop on all quadrature points
  for (auto && args : this->getArguments(el_type, ghost_type)) {
    const auto & grad_u = args["grad_u"_n];
    auto & dev_s = args["sigma_dev"_n];
    auto & h = args["history"_n];

    /// Compute the first invariant of strain
    Real Theta = grad_u.trace();

    dev_s = 2 * this->mu *
            ((grad_u + grad_u.transpose()) / 2. -
             Theta * Matrix<Real, dim, dim>::Identity() / 3.);

    h.zero();
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialStandardLinearSolidDeviatoric<dim>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  Real tau = eta / Ev;

  Real dt = this->getModel().getTimeStep();
  Real exp_dt_tau = exp(-dt / tau);
  Real exp_dt_tau_2 = exp(-.5 * dt / tau);

  /// Compute the first invariant of strain
  auto gamma_inf = E_inf / this->E;
  auto gamma_v = Ev / this->E;

  auto && arguments = this->getArguments(el_type, ghost_type);
  /// Loop on all quadrature points
  for (auto && args : arguments) {
    auto && grad_u = args["grad_u"_n];
    auto && sigma = args["sigma"_n];
    auto && dev_s = args["sigma_dev"_n];
    auto && h = args["history"_n];

    auto epsilon_d = this->template gradUToEpsilon<dim>(grad_u);
    auto Theta = epsilon_d.trace();

    epsilon_d -= Matrix<Real, dim, dim>::Identity() * Theta / 3.;

    auto U_rond_prim =
        Matrix<Real, dim, dim>::Identity() * gamma_inf * this->kpa * Theta;

    auto s = 2 * this->mu * epsilon_d;
    h = exp_dt_tau * h + exp_dt_tau_2 * (s - dev_s);
    dev_s = s;
    sigma = U_rond_prim + gamma_inf * s + gamma_v * h;
  }

  this->updateDissipatedEnergy(el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialStandardLinearSolidDeviatoric<dim>::updateDissipatedEnergy(
    ElementType el_type, GhostType ghost_type) {
  Real tau = eta / Ev;

  auto dt = this->getModel().getTimeStep();

  auto gamma_v = Ev / this->E;
  auto alpha = 1. / (2. * this->mu * gamma_v);

  for (auto && [args, dis_energy] :
       zip(this->getArguments(el_type, ghost_type),
           dissipated_energy(el_type, ghost_type))) {
    const auto & grad_u = args["grad_u"_n];
    auto & dev_s = args["sigma_dev"_n];
    auto & h = args["history"_n];

    /// Compute the first invariant of strain
    auto epsilon_d = Material::gradUToEpsilon<dim>(grad_u);

    auto Theta = epsilon_d.trace();
    epsilon_d -= Matrix<Real, dim, dim>::Identity() * Theta / 3.;

    auto q = (dev_s - h) * gamma_v;
    auto q_rate = (dev_s * gamma_v - q) / tau;

    dis_energy += ((epsilon_d - alpha * q) * q_rate * dt).sum();
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialStandardLinearSolidDeviatoric<dim>::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();
  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto && type : this->getElementFilter().elementTypes(dim, _not_ghost)) {
    de += fem.integrate(dissipated_energy(type, _not_ghost), type, _not_ghost,
                        this->getElementFilter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialStandardLinearSolidDeviatoric<dim>::getDissipatedEnergy(
    const Element & element) const {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(element.type);
  auto it = make_view(dissipated_energy(element.type, _not_ghost),
                      nb_quadrature_points)
                .begin();

  AKANTU_DEBUG_OUT();
  return fem.integrate(it[element.element], element);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialStandardLinearSolidDeviatoric<dim>::getEnergy(
    const std::string & type) {
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }
  if (type == "dissipated_sls_deviatoric") {
    return getDissipatedEnergy();
  }
  return MaterialElastic<dim>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialStandardLinearSolidDeviatoric<dim>::getEnergy(
    const std::string & energy_id, const Element & element) {
  if (energy_id == "dissipated") {
    return getDissipatedEnergy(element);
  }
  if (energy_id == "dissipated_sls_deviatoric") {
    return getDissipatedEnergy(element);
  }
  return Parent::getEnergy(energy_id, element);
}

/* -------------------------------------------------------------------------- */
template class MaterialStandardLinearSolidDeviatoric<1>;
template class MaterialStandardLinearSolidDeviatoric<2>;
template class MaterialStandardLinearSolidDeviatoric<3>;
const bool material_is_allocated_sls_deviatoric [[maybe_unused]] =
    instantiateMaterial<MaterialStandardLinearSolidDeviatoric>(
        "sls_deviatoric");

} // namespace akantu
