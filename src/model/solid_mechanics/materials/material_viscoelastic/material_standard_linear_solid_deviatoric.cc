/**
 * @file   material_standard_linear_solid_deviatoric.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Vladislav Yastrebov <vladislav.yastrebov@epfl.ch>
 *
 * @date creation: Wed May 04 2011
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Material Visco-elastic
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
#include "material_standard_linear_solid_deviatoric.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialStandardLinearSolidDeviatoric<
    dim>::MaterialStandardLinearSolidDeviatoric(SolidMechanicsModel & model,
                                                const ID & id)
    : MaterialElastic<dim>(model, id), stress_dev("stress_dev", *this),
      history_integral("history_integral", *this),
      dissipated_energy("dissipated_energy", *this) {

  AKANTU_DEBUG_IN();

  this->registerParam("Eta", eta, Real(1.), _pat_parsable | _pat_modifiable,
                      "Viscosity");
  this->registerParam("Ev", Ev, Real(1.), _pat_parsable | _pat_modifiable,
                      "Stiffness of the viscous element");
  this->registerParam("Einf", E_inf, Real(1.), _pat_readable,
                      "Stiffness of the elastic element");

  auto stress_size = dim * dim;

  this->stress_dev.initialize(stress_size);
  this->history_integral.initialize(stress_size);
  this->dissipated_energy.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialStandardLinearSolidDeviatoric<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  updateInternalParameters();
  MaterialElastic<dim>::initMaterial();

  AKANTU_DEBUG_OUT();
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
    const auto & grad_u = tuple::get<"grad_u"_h>(args);
    auto & dev_s = tuple::get<"sigma_dev"_h>(args);
    auto & h = tuple::get<"history"_h>(args);

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

  Real dt = this->model.getTimeStep();
  Real exp_dt_tau = exp(-dt / tau);
  Real exp_dt_tau_2 = exp(-.5 * dt / tau);

  Matrix<Real, dim, dim> s;
  Matrix<Real, dim, dim> epsilon_d;
  Matrix<Real, dim, dim> U_rond_prim;

  /// Compute the first invariant of strain
  auto gamma_inf = E_inf / this->E;
  auto gamma_v = Ev / this->E;

  auto && arguments = this->getArguments(el_type, ghost_type);
  /// Loop on all quadrature points
  for (auto && args : arguments) {
    auto && grad_u = tuple::get<"grad_u"_h>(args);
    auto && sigma = tuple::get<"sigma"_h>(args);
    auto && dev_s = tuple::get<"sigma_dev"_h>(args);
    auto && h = tuple::get<"history"_h>(args);

    s.zero();
    sigma.zero();

    epsilon_d = this->template gradUToEpsilon<dim>(grad_u);
    auto Theta = epsilon_d.trace();

    epsilon_d -= Matrix<Real, dim, dim>::Identity() * Theta / 3.;

    U_rond_prim =
        Matrix<Real, dim, dim>::Identity() * gamma_inf * this->kpa * Theta;

    s = 2 * this->mu * epsilon_d;
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

  Matrix<Real, dim, dim> q;
  Matrix<Real, dim, dim> q_rate;
  Matrix<Real, dim, dim> epsilon_d;

  auto dt = this->model.getTimeStep();

  auto gamma_v = Ev / this->E;
  auto alpha = 1. / (2. * this->mu * gamma_v);

  for (auto && data : zip(this->getArguments(el_type, ghost_type),
                          dissipated_energy(el_type, ghost_type))) {
    auto && args = std::get<0>(data);
    auto & dis_energy = std::get<1>(data);
    const auto & grad_u = tuple::get<"grad_u"_h>(args);
    auto & dev_s = tuple::get<"sigma_dev"_h>(args);
    auto & h = tuple::get<"history"_h>(args);

    /// Compute the first invariant of strain
    epsilon_d = Material::gradUToEpsilon<dim>(grad_u);

    auto Theta = epsilon_d.trace();
    epsilon_d -= Matrix<Real, dim, dim>::Identity() * Theta / 3.;

    q = (dev_s - h) * gamma_v;
    q_rate = (dev_s * gamma_v - q) / tau;

    dis_energy += ((epsilon_d - alpha * q) * q_rate * dt).sum();
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialStandardLinearSolidDeviatoric<dim>::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (const auto & type : this->element_filter.elementTypes(dim, _not_ghost)) {
    de +=
        this->fem.integrate(dissipated_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialStandardLinearSolidDeviatoric<dim>::getDissipatedEnergy(
    const Element & element) const {
  AKANTU_DEBUG_IN();

  auto nb_quadrature_points = this->fem.getNbIntegrationPoints(element.type);
  auto it = this->dissipated_energy(element.type, _not_ghost)
                .begin(nb_quadrature_points);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[element.element], element);
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

INSTANTIATE_MATERIAL(sls_deviatoric, MaterialStandardLinearSolidDeviatoric);

} // namespace akantu
