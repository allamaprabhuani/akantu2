/**
 * @file   material_cohesive_linear_friction_coulomb.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Mathias Lebihain <mathias.lebihain@epfl.ch>
 *
 * @date creation: Tue Jan 12 2016
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear_friction_coulomb.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialCohesiveLinearFrictionCoulomb<
    dim>::MaterialCohesiveLinearFrictionCoulomb(SolidMechanicsModel & model,
                                                const ID & id)
    : MaterialParent(model, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/** \todo Frictional contributions might be computed from the contact_tractions,
 * and delta_max at the previous time step
 */
template <UInt dim>
void MaterialCohesiveLinearFrictionCoulomb<dim>::computeTraction(
    const Array<Real> & normal, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Vector<Real> normal_opening(dim);
  Vector<Real> tangential_opening(dim);

  /// loop on each quadrature point
  for (auto && data :
       zip(make_view(this->insertion_stress(el_type, ghost_type), dim),
           make_view(this->insertion_compression(el_type, ghost_type), dim),
           this->sigma_c_eff(el_type, ghost_type),
           this->delta_c_eff(el_type, ghost_type),
           this->mu_eff(el_type, ghost_type), make_view(normal, dim),
           this->delta_max(el_type, ghost_type),
           this->damage(el_type, ghost_type),
           make_view(this->tractions(el_type, ghost_type), dim),
           make_view(this->contact_tractions(el_type, ghost_type), dim),
           make_view(this->friction_force(el_type, ghost_type), dim),
           make_view(this->opening(el_type, ghost_type), dim),
           make_view(this->contact_opening(el_type, ghost_type), dim),
           this->residual_sliding(el_type, ghost_type),
           this->residual_sliding.previous(el_type, ghost_type))) {
    auto && insertion_traction = std::get<0>(data);
    auto && insertion_compression = std::get<1>(data);
    auto && sigma_c = std::get<2>(data);
    auto && delta_c = std::get<3>(data);
    auto && mu = std::get<4>(data);
    auto && normal = std::get<5>(data);
    auto && delta_max = std::get<6>(data);
    auto && damage = std::get<7>(data);
    auto && traction = std::get<8>(data);
    auto && contact_traction = std::get<9>(data);
    auto && friction_force = std::get<10>(data);
    auto && opening = std::get<11>(data);
    auto && contact_opening = std::get<12>(data);
    auto && res_sliding = std::get<13>(data);
    auto && res_sliding_prev = std::get<14>(data);

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;

    this->computeTractionOnQuad(
        traction, opening, normal, delta_max, delta_c, insertion_traction,
        insertion_compression, sigma_c, normal_opening, tangential_opening,
        normal_opening_norm, tangential_opening_norm, damage, penetration,
        contact_traction, contact_opening);

    if (penetration) {
      /// the definition of tau_max refers to
      /// the current contact_traction computed in computeTractionOnQuad
      Real tau_max = mu * contact_traction.norm();

      // elastic sliding
      Real delta_sliding_norm =
          std::abs(tangential_opening_norm - res_sliding_prev);

      /// tau is the norm of the friction force, acting tangentially to the
      /// surface
      Real tau = std::min(this->friction_penalty * delta_sliding_norm, tau_max);

      if ((tangential_opening_norm - res_sliding_prev) < 0.0)
        tau = -tau;

      /// from tau get the x and y components of friction, to be added in the
      /// force vector
      auto tangent_unit_vector = tangential_opening / tangential_opening_norm;
      friction_force = tau * tangent_unit_vector;

      /// update residual_sliding
      res_sliding =
          tangential_opening_norm - (std::abs(tau) / this->friction_penalty);
    } else {
      friction_force.clear();
    }
    traction += friction_force;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/** \todo The frictional contribution to the tangent stiffness matrix might be
 * missing from implicit computations
 */
// template <UInt dim>
// void MaterialCohesiveLinearFrictionCoulomb<dim>::computeTangentTraction(
//     const ElementType & el_type, Array<Real> & tangent_matrix,
//     __attribute__((unused)) const Array<Real> & normal, GhostType ghost_type)
//     {
//   AKANTU_DEBUG_IN();
//
//   /// define iterators
//   auto tangent_it = tangent_matrix.begin(dim, dim);
//   auto tangent_end = tangent_matrix.end(dim, dim);
//
//   auto normal_it = this->normal.begin(dim);
//
//   auto opening_it = this->opening(el_type, ghost_type).begin(dim);
//
//   /**
//    * NB: delta_max_it points on delta_max_previous, i.e. the
//    * delta_max related to the solution of the previous incremental
//    * step
//    */
//   auto delta_max_it = this->delta_max.previous(el_type, ghost_type).begin();
//   auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
//   auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
//   auto damage_it = this->damage(el_type, ghost_type).begin();
//   auto contact_traction_it =
//       this->contact_tractions(el_type, ghost_type).begin(dim);
//   auto contact_opening_it =
//       this->contact_opening(el_type, ghost_type).begin(dim);
//   auto res_sliding_prev_it =
//       this->residual_sliding.previous(el_type, ghost_type).begin();
//
//   Vector<Real> normal_opening(dim);
//   Vector<Real> tangential_opening(dim);
//
//   for (; tangent_it != tangent_end;
//        ++tangent_it, ++normal_it, ++opening_it, ++delta_max_it, ++sigma_c_it,
//        ++delta_c_it, ++damage_it, ++contact_opening_it,
//        ++res_sliding_prev_it,
//        ++contact_traction_it) {
//
//     Real normal_opening_norm, tangential_opening_norm;
//     bool penetration;
//
//     this->computeTangentTractionOnQuad(
//         *tangent_it, *delta_max_it, *delta_c_it, *sigma_c_it, *opening_it,
//         *normal_it, normal_opening, tangential_opening, normal_opening_norm,
//         tangential_opening_norm, *damage_it, penetration,
//         *contact_opening_it);
//
//     if (penetration) {
//       // Real damage = std::min(*delta_max_it / *delta_c_it, Real(1.));
//       // Real mu = mu_max * damage;
//       Real mu = mu_max;
//
//       // Real normal_opening_norm =
//       //     std::min(contact_opening_it->dot(*normal_it), Real(0.));
//       // Real tau_max = mu * this->penalty * (std::abs(normal_opening_norm));
//
//       Real tau_max = mu * (*contact_traction_it).norm();
//
//       // elastic sliding
//       Real delta_sliding_norm =
//           std::abs(tangential_opening_norm - *res_sliding_prev_it);
//
//       /// tau is the norm of the friction force, acting tangentially to the
//       /// surface
//       Real tau = std::min(this->friction_penalty * delta_sliding_norm,
//       tau_max);
//
//       if (tau < tau_max && tau_max > Math::getTolerance()) {
//         Matrix<Real> I(dim, dim);
//         I.eye(1.);
//
//         Matrix<Real> n_outer_n(dim, dim);
//         n_outer_n.outerProduct(*normal_it, *normal_it);
//
//         Matrix<Real> nn(n_outer_n);
//         I -= nn;
//         *tangent_it += I * this->friction_penalty;
//       }
//     }
//
//     // check if the tangential stiffness matrix is symmetric
//     //    for (UInt h = 0; h < dim; ++h){
//     //      for (UInt l = h; l < dim; ++l){
//     //        if (l > h){
//     //          Real k_ls = (*tangent_it)[dim*h+l];
//     //          Real k_us =  (*tangent_it)[dim*l+h];
//     //          //          std::cout << "k_ls = " << k_ls << std::endl;
//     //          //          std::cout << "k_us = " << k_us << std::endl;
//     //          if (std::abs(k_ls) > 1e-13 && std::abs(k_us) > 1e-13){
//     //            Real error = std::abs((k_ls - k_us) / k_us);
//     //            if (error > 1e-10){
//     //	      std::cout << "non symmetric cohesive matrix" << std::endl;
//     //	      //  std::cout << "error " << error << std::endl;
//     //            }
//     //          }
//     //        }
//     //      }
//     //    }
//   }
//
//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(cohesive_linear_friction_coulomb,
                     MaterialCohesiveLinearFrictionCoulomb);

} // namespace akantu
