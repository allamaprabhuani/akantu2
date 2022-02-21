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
template <UInt spatial_dimension>
MaterialViscoelasticMaxwell<spatial_dimension>::MaterialViscoelasticMaxwell(
    SolidMechanicsModel & model, const ID & id)
    : MaterialElastic<spatial_dimension>(model, id),
      C(voigt_h::size, voigt_h::size), S(voigt_h::size, voigt_h::size),
      sigma_v("sigma_v", *this), epsilon_v("epsilon_v", *this),
      dissipated_energy("dissipated_energy", *this),
      integral("integral", *this), mechanical_work("mechanical_work", *this) {

  AKANTU_DEBUG_IN();

  this->registerParam("Einf", Einf, Real(1.), _pat_parsable | _pat_modifiable,
                      "Stiffness of the elastic element");
  this->registerParam("Eta", Eta, _pat_parsable | _pat_modifiable,
                      "Viscosity of a Maxwell element");
  this->registerParam("Ev", Ev, _pat_parsable | _pat_modifiable,
                      "Stiffness of a Maxwell element");
  this->use_previous_stress = true;
  this->use_previous_gradu = true;

  // this->dissipated_energy.initialize(1);
  // this->integral.initialize(1);
  // this->mechanical_work.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->E = Einf + Ev.norm<L_1>();
  MaterialElastic<spatial_dimension>::initMaterial();
  AKANTU_DEBUG_ASSERT(this->Eta.size() == this->Ev.size(),
                      "Eta and Ev have different dimensions! Please correct.");
  AKANTU_DEBUG_ASSERT(
      !this->finite_deformation,
      "Current material works only in infinitesimal deformations.");

  UInt stress_size = spatial_dimension * spatial_dimension;
  this->sigma_v.initialize(stress_size * this->Ev.size());
  this->sigma_v.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<
    spatial_dimension>::updateInternalParameters() {
  MaterialElastic<spatial_dimension>::updateInternalParameters();

  Real pre_mult = 1 / (1 + this->nu) / (1 - 2 * this->nu);
  UInt n = voigt_h::size;
  Real Miiii = pre_mult * (1 - this->nu);
  Real Miijj = pre_mult * this->nu;
  Real Mijij = pre_mult * 0.5 * (1 - 2 * this->nu);

  Real Siiii = 1;
  Real Siijj = -this->nu;
  Real Sijij = (2 + 2 * this->nu);

  if (spatial_dimension == 1) {
    C(0, 0) = 1;
    S(0, 0) = 1;
  } else {
    C(0, 0) = Miiii;
    S(0, 0) = Siiii;
  }
  if (spatial_dimension >= 2) {
    C(1, 1) = Miiii;
    C(0, 1) = Miijj;
    C(1, 0) = Miijj;
    C(n - 1, n - 1) = Mijij;

    S(1, 1) = Siiii;
    S(0, 1) = Siijj;
    S(1, 0) = Siijj;
    S(n - 1, n - 1) = Sijij;
  }

  if (spatial_dimension == 3) {
    C(2, 2) = Miiii;
    C(0, 2) = Miijj;
    C(1, 2) = Miijj;
    C(2, 0) = Miijj;
    C(2, 1) = Miijj;
    C(3, 3) = Mijij;
    C(4, 4) = Mijij;

    S(2, 2) = Siiii;
    S(0, 2) = Siijj;
    S(1, 2) = Siijj;
    S(2, 0) = Siijj;
    S(2, 1) = Siijj;
    S(3, 3) = Sijij;
    S(4, 4) = Sijij;
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

  Real Siiii = 1;
  Real Siijj = -this->nu;
  Real Sijij = (2 + 2 * this->nu);
  Real Sjjjj = 1 - this->nu;
  Real pre_mult_S = 1 + this->nu;

  C(0, 0) = Miiii;
  C(1, 1) = Miiii;
  C(0, 1) = Miijj;
  C(1, 0) = Miijj;
  C(n - 1, n - 1) = Mijij;

  if (this->plane_stress) {
    S(0, 0) = Siiii;
    S(1, 1) = Siiii;
    S(0, 1) = Siijj;
    S(1, 0) = Siijj;
    S(n - 1, n - 1) = Sijij;
  } else {
    S(0, 0) = pre_mult_S * Sjjjj;
    S(1, 1) = pre_mult_S * Sjjjj;
    S(0, 1) = pre_mult_S * Siijj;
    S(1, 0) = pre_mult_S * Siijj;
    S(n - 1, n - 1) = pre_mult_S * 2;
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // NOLINTNEXTLINE(bugprone-parent-virtual-call)
  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);

  auto sigma_th_it = this->sigma_th(el_type, ghost_type).begin();

  auto previous_gradu_it = this->gradu.previous(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

  auto sigma_v_it =
      this->sigma_v(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());

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
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeStressOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
    Matrix<Real> & sigma, Tensor3<Real> & sigma_v, Real sigma_th, Real dam) {

  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  Vector<Real> voigt_current_strain(voigt_h::size);
  Vector<Real> voigt_stress(voigt_h::size);
  Vector<Real> voigt_sigma_v(voigt_h::size);

  for (UInt I = 0; I < voigt_h::size; ++I) {
    Real voigt_factor = voigt_h::factors[I];
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    voigt_current_strain(I) = voigt_factor * (grad_u(i, j) + grad_u(j, i)) / 2.;
  }

  voigt_stress = (1 - dam) * this->Einf * this->C * voigt_current_strain;

  /// update a COPY of the sigma_v
  Tensor3<Real> sigma_v_copy(sigma_v);
  updateSigmaViscOnQuad(grad_u, previous_grad_u, sigma_v_copy, dam);

  for (UInt k = 0; k < this->Eta.size(); ++k) {

    for (UInt I = 0; I < voigt_h::size; ++I) {
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_sigma_v(I) = sigma_v_copy(i, j, k);
    }
    voigt_stress += voigt_sigma_v;
  }

  for (UInt I = 0; I < voigt_h::size; ++I) {
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    sigma(i, j) = sigma(j, i) =
        voigt_stress(I) + Math::kronecker(i, j) * sigma_th;
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::updateSigmaViscOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
    Tensor3<Real> & sigma_v, Real dam) {

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;

  Vector<Real> voigt_delta_strain(voigt_h::size);
  for (UInt I = 0; I < voigt_h::size; ++I) {
    Real voigt_factor = voigt_h::factors[I];
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    voigt_delta_strain(I) =
        voigt_factor * (grad_delta_u(i, j) + grad_delta_u(j, i)) / 2.;
  }

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Real E_ef_v, exp_dt_lambda;
    computeEffectiveModulus(k, E_ef_v, exp_dt_lambda, dam);
    Vector<Real> voigt_sigma_v(voigt_h::size);
    // Vector<Real> voigt_epsilon_v(voigt_h::size);

    for (UInt I = 0; I < voigt_h::size; ++I) {
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_sigma_v(I) = sigma_v(i, j, k);
    }

    voigt_sigma_v =
        exp_dt_lambda * voigt_sigma_v + E_ef_v * this->C * voigt_delta_strain;
    // voigt_epsilon_v = 1 / this->Ev(k) / (1 - dam) * this->S * voigt_sigma_v;

    for (UInt I = 0; I < voigt_h::size; ++I) {
      // Real voigt_factor = voigt_h::factors[I];
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      sigma_v(i, j, k) = sigma_v(j, i, k) = voigt_sigma_v(I);
      // epsilon_v(i, j, k) = epsilon_v(j, i, k) =
      //     voigt_epsilon_v(I) / voigt_factor;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  auto epot = this->potential_energy(el_type).begin();
  auto sigma_v_it = this->sigma_v(el_type).begin(
      spatial_dimension, spatial_dimension, this->Eta.size());
  auto epsilon_v_it = this->epsilon_v(el_type).begin(
      spatial_dimension, spatial_dimension, this->Eta.size());

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  this->computePotentialEnergyOnQuad(grad_u, *epot, *sigma_v_it, *epsilon_v_it);
  ++epot;
  ++sigma_v_it;
  ++epsilon_v_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::
    computePotentialEnergyOnQuad(const Matrix<Real> & grad_u, Real & epot,
                                 Tensor3<Real> & sigma_v,
                                 Tensor3<Real> & epsilon_v) {

  Real trace = grad_u.trace(); // trace = (\nabla u)_{kk}

  Matrix<Real> sigma(spatial_dimension, spatial_dimension);
  // \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  // u_{ij} + \nabla u_{ji})

  Real lambda = this->nu * this->Einf / ((1 + this->nu) * (1 - 2 * this->nu));
  Real mu = this->Einf / (2 * (1 + this->nu));

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      sigma(i, j) =
          (i == j) * lambda * trace + mu * (grad_u(i, j) + grad_u(j, i));
    }
  }

  epot = 0.5 * sigma.doubleDot(grad_u);

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Matrix<Real> stress_v = sigma_v(k);
    Matrix<Real> strain_v = epsilon_v(k);
    epot += 0.5 * stress_v.doubleDot(strain_v);
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::afterSolveStep(
    bool converged) {

  Material::afterSolveStep(converged);

  if (not converged) {
    return;
  }

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    // if (this->update_variable_flag) {
    auto previous_gradu_it = this->gradu.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);

    auto sigma_v_it =
        this->sigma_v(el_type, _not_ghost)
            .begin(spatial_dimension, spatial_dimension, this->Eta.size());

    // auto epsilon_v_it =
    //     this->epsilon_v(el_type, _not_ghost)
    //         .begin(spatial_dimension, spatial_dimension, this->Eta.size());

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

    updateSigmaViscOnQuad(grad_u, *previous_gradu_it, *sigma_v_it);

    ++previous_gradu_it;
    ++sigma_v_it;
    // ++epsilon_v_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
    // }
    // this->updateDissipatedEnergy(el_type);
  }
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  this->computeTangentModuliOnQuad(tangent);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, Real dam) {

  Real E_ef = this->Einf * (1 - dam);

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Real E_ef_v;
    Real exp_dt_lambda;
    this->computeEffectiveModulus(k, E_ef_v, exp_dt_lambda, dam);
    E_ef += E_ef_v;
  }

  tangent.copy(this->C);
  tangent *= E_ef;
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::updateDissipatedEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  this->computePotentialEnergy(el_type);

  auto epot = this->potential_energy(el_type).begin();
  auto dis_energy = this->dissipated_energy(el_type).begin();
  auto integral = this->integral(el_type).begin();
  auto mech_work = this->mechanical_work(el_type).begin();
  auto sigma_v_it = this->sigma_v(el_type).begin(
      spatial_dimension, spatial_dimension, this->Eta.size());
  auto epsilon_v_it = this->epsilon_v(el_type).begin(
      spatial_dimension, spatial_dimension, this->Eta.size());
  auto previous_gradu_it =
      this->gradu.previous(el_type).begin(spatial_dimension, spatial_dimension);
  auto previous_sigma_it = this->stress.previous(el_type).begin(
      spatial_dimension, spatial_dimension);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  updateDissipatedEnergyOnQuad(grad_u, *previous_gradu_it, sigma,
                               *previous_sigma_it, *dis_energy, *integral,
                               *mech_work, *epot);
  ++previous_gradu_it;
  ++previous_sigma_it;
  ++dis_energy;
  ++integral;
  ++mech_work;
  ++epot;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::
    updateDissipatedEnergyOnQuad(const Matrix<Real> & grad_u,
                                 const Matrix<Real> & previous_grad_u,
                                 const Matrix<Real> & sigma,
                                 const Matrix<Real> & previous_sigma,
                                 Real & dis_energy, Real & integral,
                                 Real & mech_work, Real pot_energy) {

  Real dt = this->model.getTimeStep();

  Matrix<Real> strain_rate = grad_u;
  strain_rate -= previous_grad_u;
  strain_rate /= dt;

  Matrix<Real> av_stress = sigma;
  av_stress += previous_sigma;
  av_stress /= 2;

  integral += av_stress.doubleDot(strain_rate) * dt;
  mech_work += std::abs(av_stress.doubleDot(strain_rate)) * dt;

  dis_energy = integral - pot_energy;
  if (std::abs(dis_energy) < std::numeric_limits<Real>::epsilon())
    dis_energy = 0;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getDissipatedEnergy()
    const {
  AKANTU_DEBUG_IN();

  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    de +=
        this->fem.integrate(this->dissipated_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getDissipatedEnergy(
    ElementType type, UInt index) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);
  auto it =
      this->dissipated_energy(type, _not_ghost).begin(nb_quadrature_points);
  UInt gindex = (this->element_filter(type, _not_ghost))(index);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[index], type, gindex);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getMechanicalWork() const {
  AKANTU_DEBUG_IN();

  Real mw = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    mw +=
        this->fem.integrate(this->mechanical_work(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return mw;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getMechanicalWork(
    ElementType type, UInt index) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);
  auto it = this->mechanical_work(type, _not_ghost).begin(nb_quadrature_points);
  UInt gindex = (this->element_filter(type, _not_ghost))(index);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[index], type, gindex);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getPotentialEnergy()
    const {
  AKANTU_DEBUG_IN();

  Real epot = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    epot +=
        this->fem.integrate(this->potential_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getPotentialEnergy(
    ElementType type, UInt index) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);
  auto it =
      this->potential_energy(type, _not_ghost).begin(nb_quadrature_points);
  UInt gindex = (this->element_filter(type, _not_ghost))(index);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[index], type, gindex);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getEnergy(
    const std::string & type) {
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }
  if (type == "potential") {
    return getPotentialEnergy();
  }
  if (type == "work") {
    return getMechanicalWork();
  }
  return MaterialElastic<spatial_dimension>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getEnergy(
    const std::string & energy_id, ElementType type, UInt index) {
  if (energy_id == "dissipated") {
    return getDissipatedEnergy(type, index);
  }
  if (energy_id == "potential") {
    return getPotentialEnergy(type, index);
  }
  if (energy_id == "work") {
    return getMechanicalWork(type, index);
  }
  return MaterialElastic<spatial_dimension>::getEnergy(energy_id, type, index);
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeEffectiveModulus(
    UInt & visc_nb, Real & E_ef_v, Real & exp_dt_lambda, Real dam) {
  Real dt = this->model.getTimeStep();
  Real lambda = (1 - dam) * this->Eta(visc_nb) / this->Ev(visc_nb) / (1 - dam);
  exp_dt_lambda = exp(-dt / lambda);

  if (Math::are_float_equal(exp_dt_lambda, 1)) {
    E_ef_v = (1 - dam) * this->Ev(visc_nb);
    // } else if (Math::are_float_equal(exp_dt_lambda, 0)) {
    //   E_ef_v = 0.;
  } else {
    E_ef_v = std::min((1 - exp_dt_lambda) * (1 - dam) * this->Ev(visc_nb) *
                          lambda / dt,
                      (1 - dam) * this->Ev(visc_nb));
  }
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(viscoelastic_maxwell, MaterialViscoelasticMaxwell);

} // namespace akantu
