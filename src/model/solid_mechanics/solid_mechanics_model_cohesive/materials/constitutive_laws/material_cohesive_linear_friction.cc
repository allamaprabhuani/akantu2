/**
 * @file   material_cohesive_linear_friction.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Tue Jan 12 2016
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_linear_friction.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveLinearFriction<spatial_dimension>::
    MaterialCohesiveLinearFriction(SolidMechanicsModel & model, const ID & id)
    : MaterialParent(model, id), residual_sliding("residual_sliding", *this),
      friction_force("friction_force", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("mu", mu_max, Real(0.), _pat_parsable | _pat_readable,
                      "Maximum value of the friction coefficient");

  this->registerParam("penalty_for_friction", friction_penalty, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Penalty parameter for the friction behavior");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFriction<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialParent::initMaterial();

  friction_force.initialize(spatial_dimension);
  residual_sliding.initialize(1);
  residual_sliding.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
// template <UInt spatial_dimension>
// void MaterialCohesiveLinearFriction<spatial_dimension>::computeTraction(
//     __attribute__((unused)) const Array<Real> & normal, ElementType el_type,
//     GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   residual_sliding.resize();
//   friction_force.resize();

//   auto & traction = this->tractions(el_type, ghost_type);
//   auto & opening = this->opening(el_type, ghost_type);
//   auto & opening_prev = this->opening.previous(el_type, ghost_type);
//   auto & contact_traction = this->contact_tractions(el_type, ghost_type);
//   auto & contact_opening = this->contact_opening(el_type, ghost_type);
//   auto & sigma_c = this->sigma_c_eff(el_type, ghost_type);
//   auto & delta_max = this->delta_max(el_type, ghost_type);
//   auto & delta_max_prev = this->delta_max.previous(el_type, ghost_type);
//   auto & delta_c = this->delta_c_eff(el_type, ghost_type);
//   auto & damage = this->damage(el_type, ghost_type);
//   auto & insertion_stress = this->insertion_stress(el_type, ghost_type);
//   auto & res_sliding = this->residual_sliding(el_type, ghost_type);
//   auto & res_sliding_prev =
//       this->residual_sliding.previous(el_type, ghost_type);
//   auto & friction_force = this->friction_force(el_type, ghost_type);
//   auto & penetration = this->penetration(el_type, ghost_type);

//   Vector<Real> normal_opening(spatial_dimension);
//   Vector<Real> tangential_opening(spatial_dimension);

//   if (not this->model->isDefaultSolverExplicit()) {
//     this->delta_max(el_type, ghost_type)
//         .copy(this->delta_max.previous(el_type, ghost_type));
//   }

//   /// loop on each quadrature point
//   for (auto && data :
//        zip(make_view(traction, spatial_dimension),
//            make_view(opening, spatial_dimension),
//            make_view(opening_prev, spatial_dimension),
//            make_view(normal, spatial_dimension), sigma_c, delta_max,
//            delta_max_prev, delta_c, damage,
//            make_view(insertion_stress, spatial_dimension), res_sliding,
//            res_sliding_prev, make_view(friction_force, spatial_dimension),
//            penetration, make_view(contact_opening, spatial_dimension),
//            make_view(contact_traction, spatial_dimension))) {
//     auto & _traction = std::get<0>(data);
//     auto & _opening = std::get<1>(data);
//     auto & _opening_prev = std::get<2>(data);
//     auto & _normal = std::get<3>(data);
//     auto & _sigma_c = std::get<4>(data);
//     auto & _delta_max = std::get<5>(data);
//     auto & _delta_max_prev = std::get<6>(data);
//     auto & _delta_c = std::get<7>(data);
//     auto & _damage = std::get<8>(data);
//     auto & _insertion_stress = std::get<9>(data);
//     auto & _res_sliding = std::get<10>(data);
//     auto & _res_sliding_prev = std::get<11>(data);
//     auto & _friction_force = std::get<12>(data);
//     auto & _penetration = std::get<13>(data);
//     auto & _contact_opening = std::get<14>(data);
//     auto & _contact_traction = std::get<15>(data);
//     Real normal_opening_norm;
//     Real tangential_opening_norm;

//     this->computeTractionOnQuad(
//         _traction, _opening, _normal, _delta_max, _delta_c,
//         _insertion_stress, _sigma_c, normal_opening, tangential_opening,
//         normal_opening_norm, tangential_opening_norm, _damage, _penetration,
//         _contact_traction, _contact_opening);

//     if (_penetration) {
//       /// the friction coefficient mu increases with the damage. It
//       /// equals the maximum value when damage = 1.
//       //      Real damage = std::min(*delta_max_prev_it / *delta_c_it,
//       //      Real(1.));
//       Real mu = mu_max; // * damage;
//       // Real contact_opening_norm =
//       //     std::min(_contact_opening.dot(_normal), Real(0.));
//       Real contact_opening_norm = _contact_opening.dot(_normal);

//       /// the definition of tau_max refers to the opening
//       /// (penetration) of the previous incremental step

//       Real tau_max = mu * this->penalty * (std::abs(contact_opening_norm));
//       Real trial_elastic_slip = tangential_opening_norm - _res_sliding_prev;

//       /// tau is the norm of the friction force, acting tangentially to the
//       /// surface
//       Real tau =
//           std::min(std::abs(friction_penalty * trial_elastic_slip), tau_max);

//       if (trial_elastic_slip < 0.0) {
//         tau = -tau;
//       }

//       /// from tau get the x and y components of friction, to be added in the
//       /// force vector
//       if (tangential_opening_norm > Math::getTolerance()) {
//         Vector<Real> tangent_unit_vector(spatial_dimension);
//         tangent_unit_vector = tangential_opening / tangential_opening_norm;
//         _friction_force = tau * tangent_unit_vector;
//       } else {
//         _friction_force.zero();
//       }
//       /// update residual_sliding
//       if (friction_penalty == 0.) {
//         _res_sliding = tangential_opening_norm;
//       } else {
//         _res_sliding =
//             tangential_opening_norm - (std::abs(tau) / friction_penalty);
//       }
//     } else {
//       _friction_force.zero();
//     }

//     _traction += _friction_force;
//   }

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFriction<spatial_dimension>::computeTraction(
    __attribute__((unused)) const Array<Real> & normal, ElementType el_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  residual_sliding.resize();
  friction_force.resize();

  auto & traction = this->tractions(el_type, ghost_type);
  auto & opening = this->opening(el_type, ghost_type);
  auto & opening_prev = this->opening.previous(el_type, ghost_type);
  auto & contact_traction = this->contact_tractions(el_type, ghost_type);
  auto & contact_opening = this->contact_opening(el_type, ghost_type);
  auto & sigma_c = this->sigma_c_eff(el_type, ghost_type);
  auto & delta_max = this->delta_max(el_type, ghost_type);
  auto & delta_max_prev = this->delta_max.previous(el_type, ghost_type);
  auto & delta_c = this->delta_c_eff(el_type, ghost_type);
  auto & damage = this->damage(el_type, ghost_type);
  auto & insertion_stress = this->insertion_stress(el_type, ghost_type);
  auto & res_sliding = this->residual_sliding(el_type, ghost_type);
  auto & res_sliding_prev =
      this->residual_sliding.previous(el_type, ghost_type);
  auto & friction_force = this->friction_force(el_type, ghost_type);
  auto & penetration = this->penetration(el_type, ghost_type);

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  if (not this->model->isDefaultSolverExplicit()) {
    this->delta_max(el_type, ghost_type)
        .copy(this->delta_max.previous(el_type, ghost_type));
  }

  /// loop on each quadrature point
  for (auto && data :
       zip(make_view(traction, spatial_dimension),
           make_view(opening, spatial_dimension),
           make_view(opening_prev, spatial_dimension),
           make_view(normal, spatial_dimension), sigma_c, delta_max,
           delta_max_prev, delta_c, damage,
           make_view(insertion_stress, spatial_dimension), res_sliding,
           res_sliding_prev, make_view(friction_force, spatial_dimension),
           penetration, make_view(contact_opening, spatial_dimension),
           make_view(contact_traction, spatial_dimension))) {
    auto & _traction = std::get<0>(data);
    auto & _opening = std::get<1>(data);
    auto & _opening_prev = std::get<2>(data);
    auto & _normal = std::get<3>(data);
    auto & _sigma_c = std::get<4>(data);
    auto & _delta_max = std::get<5>(data);
    auto & _delta_max_prev = std::get<6>(data);
    auto & _delta_c = std::get<7>(data);
    auto & _damage = std::get<8>(data);
    auto & _insertion_stress = std::get<9>(data);
    auto & _res_sliding = std::get<10>(data);
    auto & _res_sliding_prev = std::get<11>(data);
    auto & _friction_force = std::get<12>(data);
    auto & _penetration = std::get<13>(data);
    auto & _contact_opening = std::get<14>(data);
    auto & _contact_traction = std::get<15>(data);
    Real normal_opening_norm;
    Real tangential_opening_norm;

    this->computeTractionOnQuad(
        _traction, _opening, _normal, _delta_max, _delta_c, _insertion_stress,
        _sigma_c, normal_opening, tangential_opening, normal_opening_norm,
        tangential_opening_norm, _damage, _penetration, _contact_traction,
        _contact_opening);

    if (_penetration) {
      /// the friction coefficient mu increases with the damage. It
      /// equals the maximum value when damage = 1.
      //      Real damage = std::min(*delta_max_prev_it / *delta_c_it,
      //      Real(1.));
      Real mu = mu_max; // * damage;
      // Real contact_opening_norm =
      //     std::min(_contact_opening.dot(_normal), Real(0.));
      Real contact_opening_norm = _contact_opening.dot(_normal);

      /// the definition of tau_max refers to the opening
      /// (penetration) of the previous incremental step

      auto tau_max = mu * this->penalty * (std::abs(contact_opening_norm));
      auto trial_elastic_slip = tangential_opening_norm - _res_sliding_prev;

      auto tau_trial = std::abs(friction_penalty * trial_elastic_slip);
      auto tau = std::min(tau_trial, tau_max);

      if (trial_elastic_slip < 0.0) {
        tau = -tau;
      }

      /// from tau get the x and y components of friction, to be added in the
      /// force vector
      if (tangential_opening_norm > Math::getTolerance()) {
        Vector<Real> tangent_unit_vector(spatial_dimension);
        tangent_unit_vector = tangential_opening / tangential_opening_norm;
        _friction_force = tau * tangent_unit_vector;
      } else {
        _friction_force.zero();
      }
      /// update residual_sliding
      if (tau_trial > tau_max) {
        if (friction_penalty == 0.) {
          _res_sliding = tangential_opening_norm;
        } else {
          auto lambda = (tau_trial - tau_max) / friction_penalty;
          _res_sliding = _res_sliding_prev + lambda;
        }
      } else {
        _res_sliding = _res_sliding_prev;
      }

    } else {
      _friction_force.zero();
    }

    _traction += _friction_force;
  }

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFriction<spatial_dimension>::computeTangentTraction(
    ElementType el_type, Array<Real> & tangent_matrix,
    __attribute__((unused)) const Array<Real> & normal, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & opening = this->opening(el_type, ghost_type);
  auto & opening_prev = this->opening.previous(el_type, ghost_type);
  auto & delta_max_prev = this->delta_max.previous(el_type, ghost_type);
  auto & sigma_c = this->sigma_c_eff(el_type, ghost_type);
  auto & delta_c = this->delta_c_eff(el_type, ghost_type);
  auto & damage = this->damage(el_type, ghost_type);
  auto & contact_opening = this->contact_opening(el_type, ghost_type);
  auto & contact_opening_prev =
      this->contact_opening.previous(el_type, ghost_type);
  auto & res_sliding_prev =
      this->residual_sliding.previous(el_type, ghost_type);
  auto & penetration = this->penetration(el_type, ghost_type);

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  for (auto && data :
       zip(make_view(tangent_matrix, spatial_dimension, spatial_dimension),
           make_view(opening, spatial_dimension),
           make_view(opening_prev, spatial_dimension),
           make_view(contact_opening, spatial_dimension),
           make_view(contact_opening_prev, spatial_dimension),
           make_view(normal, spatial_dimension), sigma_c, delta_max_prev,
           delta_c, damage, res_sliding_prev, penetration)) {
    auto & _tangent = std::get<0>(data);
    auto & _opening = std::get<1>(data);
    auto & _opening_prev = std::get<2>(data);
    auto & _contact_opening = std::get<3>(data);
    auto & _contact_opening_prev = std::get<4>(data);
    auto & _normal = std::get<5>(data);
    auto & _sigma_c = std::get<6>(data);
    auto & _delta_max_prev = std::get<7>(data);
    auto & _delta_c = std::get<8>(data);
    auto & _damage = std::get<9>(data);
    auto & _res_sliding_prev = std::get<10>(data);
    auto & _penetration = std::get<11>(data);

    Real normal_opening_norm;
    Real tangential_opening_norm;
    this->computeTangentTractionOnQuad(
        _tangent, _delta_max_prev, _delta_c, _sigma_c, _opening, _normal,
        normal_opening, tangential_opening, normal_opening_norm,
        tangential_opening_norm, _damage, _penetration, _contact_opening);

    if (_penetration) {
      //      Real damage = std::min(*delta_max_it / *delta_c_it, Real(1.));
      Real mu = mu_max; // * damage;

      Real contact_opening_norm =
          std::min(_contact_opening.dot(_normal), Real(0.));
      //      Vector<Real> normal_opening_prev = (*normal_it);
      //      normal_opening_prev *= normal_opening_prev_norm;

      Real tau_max = mu * this->penalty * (std::abs(contact_opening_norm));
      Real trial_elastic_slip = tangential_opening_norm - _res_sliding_prev;

      // tau is the norm of the friction force, acting tangentially to the
      // surface
      Real tau =
          std::min(std::abs(friction_penalty * trial_elastic_slip), tau_max);

      if ((tau < tau_max * 1.05 && tau_max > Math::getTolerance()) or
          Math::are_float_equal(tau_max, 0.)) {
        Matrix<Real> I(spatial_dimension, spatial_dimension);
        I.eye(1.);

        Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
        n_outer_n.outerProduct(_normal, _normal);

        Matrix<Real> nn(n_outer_n);
        I -= nn;
        _tangent += I * friction_penalty;
      }
    }

    // check if the tangential stiffness matrix is symmetric
    //    for (UInt h = 0; h < spatial_dimension; ++h){
    //      for (UInt l = h; l < spatial_dimension; ++l){
    //        if (l > h){
    //          Real k_ls = (*tangent_it)[spatial_dimension*h+l];
    //          Real k_us =  (*tangent_it)[spatial_dimension*l+h];
    //          //          std::cout << "k_ls = " << k_ls << std::endl;
    //          //          std::cout << "k_us = " << k_us << std::endl;
    //          if (std::abs(k_ls) > 1e-13 && std::abs(k_us) > 1e-13){
    //            Real error = std::abs((k_ls - k_us) / k_us);
    //            if (error > 1e-10){
    //	      std::cout << "non symmetric cohesive matrix" << std::endl;
    //	      //  std::cout << "error " << error << std::endl;
    //            }
    //          }
    //        }
    //      }
    //    }
  }

  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */

INSTANTIATE_MATERIAL(cohesive_linear_friction, MaterialCohesiveLinearFriction);

} // namespace akantu
