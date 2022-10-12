/**
 * @file   material_viscoelastic_maxwell.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Mon Jun 04 2018
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Material Visco-elastic, based on Maxwell chain,
 * see
 * [] R. de Borst and A.H. van den Boogaard "Finite-element modeling of
 * deformation and cracking in early-age concrete", J.Eng.Mech., 1994
 * as well as
 * [] Manual of DIANA FEA Theory manual v.10.2 Section 37.6
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_viscoelastic_maxwell.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialViscoelasticMaxwell<dim>::MaterialViscoelasticMaxwell(
    SolidMechanicsModel & model, const ID & id)
    : MaterialElastic<dim>(model, id), C(voigt_h::size, voigt_h::size),
      D(voigt_h::size, voigt_h::size), sigma_v("sigma_v", *this),
      epsilon_v("epsilon_v", *this),
      dissipated_energy("dissipated_energy", *this),
      mechanical_work("mechanical_work", *this) {

  AKANTU_DEBUG_IN();

  this->registerParam("Einf", Einf, Real(1.), _pat_parsable | _pat_modifiable,
                      "Stiffness of the elastic element");
  this->registerParam("previous_dt", previous_dt, Real(0.), _pat_readable,
                      "Time step of previous solveStep");
  this->registerParam("Eta", Eta, _pat_parsable | _pat_modifiable,
                      "Viscosity of a Maxwell element");
  this->registerParam("Ev", Ev, _pat_parsable | _pat_modifiable,
                      "Stiffness of a Maxwell element");
  this->update_variable_flag = true;
  this->use_previous_stress = true;
  this->use_previous_gradu = true;

  this->dissipated_energy.initialize(1);
  this->mechanical_work.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialViscoelasticMaxwell<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->E = Einf + Ev.lpNorm<1>();
  //  this->E = std::min(this->Einf, this->Ev(0));
  MaterialElastic<dim>::initMaterial();

  AKANTU_DEBUG_ASSERT(this->Eta.size() == this->Ev.size(),
                      "Eta and Ev have different dimensions! Please correct.");
  AKANTU_DEBUG_ASSERT(
      !this->finite_deformation,
      "Current material works only in infinitesimal deformations.");

  UInt stress_size = dim * dim;
  this->sigma_v.initialize(stress_size * this->Ev.size());
  this->epsilon_v.initialize(stress_size * this->Ev.size());

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::updateInternalParameters() {
  MaterialElastic<dim>::updateInternalParameters();

  Real pre_mult = 1 / (1 + this->nu) / (1 - 2 * this->nu);
  UInt n = voigt_h::size;
  Real Miiii = pre_mult * (1 - this->nu);
  Real Miijj = pre_mult * this->nu;
  Real Mijij = pre_mult * 0.5 * (1 - 2 * this->nu);

  Real Diiii = 1;
  Real Diijj = -this->nu;
  Real Dijij = (2 + 2 * this->nu);

  if (dim == 1) {
    C(0, 0) = 1;
    D(0, 0) = 1;
  } else {
    C(0, 0) = Miiii;
    D(0, 0) = Diiii;
  }
  if (dim >= 2) {
    C(1, 1) = Miiii;
    C(0, 1) = Miijj;
    C(1, 0) = Miijj;
    C(n - 1, n - 1) = Mijij;

    D(1, 1) = Diiii;
    D(0, 1) = Diijj;
    D(1, 0) = Diijj;
    D(n - 1, n - 1) = Dijij;
  }

  if (dim == 3) {
    C(2, 2) = Miiii;
    C(0, 2) = Miijj;
    C(1, 2) = Miijj;
    C(2, 0) = Miijj;
    C(2, 1) = Miijj;
    C(3, 3) = Mijij;
    C(4, 4) = Mijij;

    D(2, 2) = Diiii;
    D(0, 2) = Diijj;
    D(1, 2) = Diijj;
    D(2, 0) = Diijj;
    D(2, 1) = Diijj;
    D(3, 3) = Dijij;
    D(4, 4) = Dijij;
  }
}

/* -------------------------------------------------------------------------- */
template <> void MaterialViscoelasticMaxwell<2>::updateInternalParameters() {
  MaterialElastic<2>::updateInternalParameters();

  Real pre_mult = 1 / (1 + this->nu) / (1 - 2 * this->nu);
  UInt n = voigt_h::size;
  Real Miiii = pre_mult * (1 - this->nu);
  Real Miijj = pre_mult * this->nu;
  Real Mijij = pre_mult * 0.5 * (1 - 2 * this->nu);

  Real Diiii = 1;
  Real Diijj = -this->nu;
  Real Dijij = (2 + 2 * this->nu);

  C(0, 0) = Miiii;
  C(1, 1) = Miiii;
  C(0, 1) = Miijj;
  C(1, 0) = Miijj;
  C(n - 1, n - 1) = Mijij;

  D(0, 0) = Diiii;
  D(1, 1) = Diiii;
  D(0, 1) = Diijj;
  D(1, 0) = Diijj;
  D(n - 1, n - 1) = Dijij;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::computeStress(ElementType el_type,
                                                     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // NOLINTNEXTLINE(bugprone-parent-virtual-call)
  MaterialThermal<dim>::computeStress(el_type, ghost_type);

  auto sigma_th_it = this->sigma_th(el_type, ghost_type).begin();

  auto previous_gradu_it =
      make_view<dim, dim>(this->gradu.previous(el_type, ghost_type)).begin();

  auto previous_stress_it =
      make_view<dim, dim>(this->stress.previous(el_type, ghost_type)).begin();

  auto sigma_v_it =
      make_view(this->sigma_v(el_type, ghost_type), dim, dim, this->Eta.size())
          .begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeStressOnQuad(grad_u, *previous_gradu_it, sigma, *sigma_v_it,
                      *sigma_th_it);

  ++sigma_th_it;
  ++previous_gradu_it;
  ++sigma_v_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1, typename D2, typename D3>
void MaterialViscoelasticMaxwell<dim>::computeStressOnQuad(
    const Eigen::MatrixBase<D1> & grad_u,
    const Eigen::MatrixBase<D2> & previous_grad_u,
    Eigen::MatrixBase<D3> & sigma, Tensor3Proxy<Real> & sigma_v,
    const Real & sigma_th) {
  Real dt = this->model.getTimeStep();
  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  auto voigt_current_strain =
      Material::strainToVoigt<dim>(Material::gradUToEpsilon<dim>(grad_u));
  auto voigt_previous_strain = Material::strainToVoigt<dim>(
      Material::gradUToEpsilon<dim>(previous_grad_u));

  Vector<Real, voigt_h::size> voigt_stress =
      this->Einf * this->C * voigt_current_strain;

  Vector<Real, voigt_h::size> stress =
      this->C * (voigt_current_strain - voigt_previous_strain);

  for (auto && [Eta, Ev, sigma] : zip(this->Eta, this->Ev, sigma_v)) {
    auto lambda = Eta / Ev;
    auto exp_dt_lambda = exp(-dt / lambda);
    Real E_additional;

    if (exp_dt_lambda == 1) {
      E_additional = Ev;
    } else {
      E_additional = (1 - exp_dt_lambda) * Ev * lambda / dt;
    }

    voigt_stress += E_additional * stress +
                    exp_dt_lambda * Material::stressToVoigt<dim>(sigma);
  }

  for (Int I = 0; I < voigt_h::size; ++I) {
    auto i = voigt_h::vec[I][0];
    auto j = voigt_h::vec[I][1];

    sigma(i, j) = sigma(j, i) =
        voigt_stress(I) + Math::kronecker(i, j) * sigma_th;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  auto epot = this->potential_energy(el_type).begin();
  auto sigma_v_it = this->sigma_v(el_type).begin(dim, dim, this->Eta.size());
  auto epsilon_v_it =
      this->epsilon_v(el_type).begin(dim, dim, this->Eta.size());

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  this->computePotentialEnergyOnQuad(grad_u, *epot, *sigma_v_it, *epsilon_v_it);

  ++epot;
  ++sigma_v_it;
  ++epsilon_v_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1>
void MaterialViscoelasticMaxwell<dim>::computePotentialEnergyOnQuad(
    const Eigen::MatrixBase<D1> & grad_u, Real & epot,
    Tensor3Proxy<Real> & sigma_v, Tensor3Proxy<Real> & epsilon_v) {

  auto voigt_strain =
      Material::strainToVoigt<dim>(Material::gradUToEpsilon<dim>(grad_u));
  Vector<Real, voigt_h::size> voigt_stress =
      this->Einf * this->C * voigt_strain;

  epot = 0.5 * voigt_stress.dot(voigt_strain);

  for (Int k = 0; k < this->Eta.size(); ++k) {
    Matrix<Real> stress_v = sigma_v(k);
    Matrix<Real> strain_v = epsilon_v(k);
    epot += 0.5 * stress_v.doubleDot(strain_v);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::afterSolveStep(bool converged) {
  Material::afterSolveStep(converged);

  if (not converged) {
    return;
  }

  for (const auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    if (this->update_variable_flag) {
      auto previous_gradu_it =
          this->gradu.previous(el_type, _not_ghost).begin(dim, dim);

      auto sigma_v_it =
          this->sigma_v(el_type, _not_ghost).begin(dim, dim, this->Eta.size());

      auto epsilon_v_it = this->epsilon_v(el_type, _not_ghost)
                              .begin(dim, dim, this->Eta.size());

      MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

      updateIntVarOnQuad(grad_u, *previous_gradu_it, *sigma_v_it,
                         *epsilon_v_it);

      ++previous_gradu_it;
      ++sigma_v_it;
      ++epsilon_v_it;

      MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
    }
    this->updateDissipatedEnergy(el_type);
  }
}
/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1, typename D2>
void MaterialViscoelasticMaxwell<dim>::updateIntVarOnQuad(
    const Eigen::MatrixBase<D1> & grad_u,
    const Eigen::MatrixBase<D2> & previous_grad_u, Tensor3Proxy<Real> & sigma_v,
    Tensor3Proxy<Real> & epsilon_v) {

  Matrix<Real, dim, dim> grad_delta_u = grad_u - previous_grad_u;

  Real dt = this->model.getTimeStep();

  auto voigt_delta_strain =
      Material::strainToVoigt<dim>(Material::gradUToEpsilon<dim>(grad_delta_u));

  for (Idx k = 0; k < this->Eta.size(); ++k) {
    auto lambda = this->Eta(k) / this->Ev(k);
    auto exp_dt_lambda = exp(-dt / lambda);
    Real E_ef_v;

    if (exp_dt_lambda == 1) {
      E_ef_v = this->Ev(k);
    } else {
      E_ef_v = (1 - exp_dt_lambda) * this->Ev(k) * lambda / dt;
    }

    auto voigt_sigma_v = Material::stressToVoigt<dim>(sigma_v(k));

    Vector<Real, voigt_h::size> voigt_epsilon_v =
        exp_dt_lambda * voigt_sigma_v + E_ef_v * this->C * voigt_delta_strain;
    voigt_epsilon_v = 1 / Ev(k) * this->D * voigt_sigma_v;

    for (Int I = 0; I < voigt_h::size; ++I) {
      auto i = voigt_h::vec[I][0];
      auto j = voigt_h::vec[I][1];

      sigma_v(i, j, k) = sigma_v(j, i, k) = voigt_sigma_v(I);
      epsilon_v(i, j, k) = epsilon_v(j, i, k) = voigt_epsilon_v(I);
    }
  }
}
/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real dt = this->model.getTimeStep();
  Real E_ef = this->Einf;

  for (Int k = 0; k < Eta.size(); ++k) {
    Real lambda = this->Eta(k) / this->Ev(k);
    Real exp_dt_lambda = exp(-dt / lambda);
    if (exp_dt_lambda == 1) {
      E_ef += this->Ev(k);
    } else {
      E_ef += (1 - exp_dt_lambda) * this->Ev(k) * lambda / dt;
    }
  }

  this->previous_dt = dt;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  this->computeTangentModuliOnQuad(tangent);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  tangent_matrix *= E_ef;

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1>
void MaterialViscoelasticMaxwell<dim>::computeTangentModuliOnQuad(
    Eigen::MatrixBase<D1> & tangent) {

  tangent = C;
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialViscoelasticMaxwell<dim>::updateIntVariables() {

  for (const auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {

    auto previous_gradu_it =
        this->gradu.previous(el_type, _not_ghost).begin(dim, dim);

    auto sigma_v_it =
        this->sigma_v(el_type, _not_ghost).begin(dim, dim, this->Eta.size());

    auto epsilon_v_it =
        this->epsilon_v(el_type, _not_ghost).begin(dim, dim, this->Eta.size());

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

    updateIntVarOnQuad(grad_u, *previous_gradu_it, *sigma_v_it, *epsilon_v_it);

    ++previous_gradu_it;
    ++sigma_v_it;
    ++epsilon_v_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::updateDissipatedEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  this->computePotentialEnergy(el_type);

  auto epot = this->potential_energy(el_type).begin();
  auto dis_energy = this->dissipated_energy(el_type).begin();
  auto mech_work = this->mechanical_work(el_type).begin();
  auto sigma_v_it = this->sigma_v(el_type).begin(dim, dim, this->Eta.size());
  auto epsilon_v_it =
      this->epsilon_v(el_type).begin(dim, dim, this->Eta.size());
  auto previous_gradu_it = this->gradu.previous(el_type).begin(dim, dim);
  auto previous_sigma_it = this->stress.previous(el_type).begin(dim, dim);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  updateDissipatedEnergyOnQuad(grad_u, *previous_gradu_it, sigma,
                               *previous_sigma_it, *dis_energy, *mech_work,
                               *epot);
  ++previous_gradu_it;
  ++previous_sigma_it;
  ++dis_energy;
  ++mech_work;
  ++epot;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1, typename D2, typename D3, typename D4>
void MaterialViscoelasticMaxwell<dim>::updateDissipatedEnergyOnQuad(
    const Eigen::MatrixBase<D1> & grad_u,
    const Eigen::MatrixBase<D2> & previous_grad_u,
    const Eigen::MatrixBase<D3> & sigma,
    const Eigen::MatrixBase<D4> & previous_sigma, Real & dis_energy,
    Real & mech_work, const Real & pot_energy) {

  Real dt = this->model.getTimeStep();

  Matrix<Real> strain_rate = grad_u;
  strain_rate -= previous_grad_u;
  strain_rate /= dt;

  Matrix<Real> av_stress = sigma;
  av_stress += previous_sigma;
  av_stress /= 2;

  mech_work += av_stress.doubleDot(strain_rate) * dt;

  dis_energy = mech_work - pot_energy;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (const auto & type : this->element_filter.elementTypes(dim, _not_ghost)) {
    de +=
        this->fem.integrate(this->dissipated_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getDissipatedEnergy(
    const Element & element) const {
  AKANTU_DEBUG_IN();

  auto nb_quadrature_points = this->fem.getNbIntegrationPoints(element.type);
  auto it = this->dissipated_energy(element.type, _not_ghost)
                .begin(nb_quadrature_points);
  auto gindex = this->element_filter(element);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[element.element],
                             {element.type, gindex, _not_ghost});
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getMechanicalWork() const {
  AKANTU_DEBUG_IN();

  Real mw = 0.;

  /// integrate the dissipated energy for each type of elements
  for (const auto & type : this->element_filter.elementTypes(dim, _not_ghost)) {
    mw +=
        this->fem.integrate(this->mechanical_work(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return mw;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getMechanicalWork(
    const Element & element) const {
  AKANTU_DEBUG_IN();

  auto nb_quadrature_points = this->fem.getNbIntegrationPoints(element.type);
  auto it = this->mechanical_work(element.type, _not_ghost)
                .begin(nb_quadrature_points);
  auto gindex = this->element_filter(element);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[element.element],
                             {element.type, gindex, _not_ghost});
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getPotentialEnergy() const {
  AKANTU_DEBUG_IN();

  Real epot = 0.;

  /// integrate the dissipated energy for each type of elements
  for (const auto & type : this->element_filter.elementTypes(dim, _not_ghost)) {
    epot +=
        this->fem.integrate(this->potential_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getPotentialEnergy(
    const Element & element) const {
  AKANTU_DEBUG_IN();

  auto nb_quadrature_points = this->fem.getNbIntegrationPoints(element.type);
  auto it = this->potential_energy(element.type, _not_ghost)
                .begin(nb_quadrature_points);
  auto gindex = this->element_filter(element);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[element.element],
                             {element.type, gindex, _not_ghost});
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getEnergy(const std::string & type) {
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }
  if (type == "potential") {
    return getPotentialEnergy();
  }
  if (type == "work") {
    return getMechanicalWork();
  }
  return MaterialElastic<dim>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticMaxwell<dim>::getEnergy(const std::string & energy_id,
                                                 const Element & element) {
  if (energy_id == "dissipated") {
    return getDissipatedEnergy(element);
  }
  if (energy_id == "potential") {
    return getPotentialEnergy(element);
  }
  if (energy_id == "work") {
    return getMechanicalWork(element);
  }
  return MaterialElastic<dim>::getEnergy(energy_id, element);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::forceUpdateVariable() {
  update_variable_flag = true;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticMaxwell<dim>::forceNotUpdateVariable() {
  update_variable_flag = false;
}

/* -------------------------------------------------------------------------- */
template class MaterialViscoelasticMaxwell<1>;
template class MaterialViscoelasticMaxwell<2>;
template class MaterialViscoelasticMaxwell<3>;
static bool material_is_allocated_viscoelastic_maxwell =
    instantiateMaterial<MaterialViscoelasticMaxwell>("viscoelastic_maxwell");

} // namespace akantu
