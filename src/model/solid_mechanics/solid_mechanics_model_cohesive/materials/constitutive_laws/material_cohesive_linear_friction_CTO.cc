/**
 * @file   material_cohesive_linear_friction_CTO.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@gmail.com>
 *
 * @date creation: Thu Sept 21 2023
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Addition to the existing linear friction cohesive material
 * by adding the consistent tangent operator when in sliping mode
 * Implementation is based on the P. Wriggers book on Computational Contact
 * Mechanics 2nd ed, subchapter 8.2, equation 8.45. The only difference with his
 * implementation is the fact that at the level of cohesive element we assemble
 * a tangent operator as if it would be a solid element. The part of computing
 * jumps in displacement between pairs of opposite nodes is taken in a different
 * part of Akantu. Meanwhile Wriggers includes this jump-computing operation
 * into the operator itself. That's why vectors Ci and T_alpha consist of {n1,
 * -n1} and {e_alhpa, -e_alpha} correspondignly, while for us, those are simply
 * {n1} and {e_alpha}.
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
#include "material_cohesive_linear_friction_CTO.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveLinearFrictionCTO<spatial_dimension>::
    MaterialCohesiveLinearFrictionCTO(SolidMechanicsModel & model,
                                      const ID & id)
    : MaterialParent(model, id), basis("basis", *this) {}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFrictionCTO<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialParent::initMaterial();

  basis.initialize(spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFrictionCTO<
    spatial_dimension>::computeTangentTraction(ElementType el_type,
                                               Array<Real> & tangent_matrix,
                                               __attribute__((unused))
                                               const Array<Real> & normal,
                                               GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & basis = this->basis(el_type, ghost_type);
  auto && position = this->model->getMesh().getNodes();
  this->computeBasis(position, basis, el_type, ghost_type);

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
           delta_c, damage, res_sliding_prev, penetration,
           make_view(basis, spatial_dimension, spatial_dimension))) {
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
    const auto & _basis = std::get<12>(data);
    Real normal_opening_norm;
    Real tangential_opening_norm;
    this->computeTangentTractionOnQuad(
        _tangent, _delta_max_prev, _delta_c, _sigma_c, _opening, _normal,
        normal_opening, tangential_opening, normal_opening_norm,
        tangential_opening_norm, _damage, _penetration, _contact_opening);

    if (_penetration) {
      //      Real damage = std::min(*delta_max_it / *delta_c_it, Real(1.));
      Real mu = this->mu_max; // * damage;

      Real contact_opening_norm =
          std::min(_contact_opening.dot(_normal), Real(0.));
      //      Vector<Real> normal_opening_prev = (*normal_it);
      //      normal_opening_prev *= normal_opening_prev_norm;

      Real tau_max = mu * this->penalty * (std::abs(contact_opening_norm));
      Real trial_elastic_slip = tangential_opening_norm - _res_sliding_prev;

      // tau is the norm of the friction force, acting tangentially to the
      // surface
      Real tau = std::min(std::abs(this->friction_penalty * trial_elastic_slip),
                          tau_max);

      if ((tau < tau_max && tau_max > Math::getTolerance()) or
          Math::are_float_equal(tau_max, 0.)) {
        Matrix<Real> I(spatial_dimension, spatial_dimension);
        I.eye(1.);

        Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
        n_outer_n.outerProduct(_normal, _normal);

        Matrix<Real> nn(n_outer_n);
        I -= nn;
        _tangent += I * this->friction_penalty;
      } else if (tau == tau_max && tau_max > Math::getTolerance()) {
        AKANTU_DEBUG_ASSERT(
            spatial_dimension == 3,
            "Consistent Tangent operator is implemented only in 3D");

        Vector<Real> sliding_dir(tangential_opening);
        sliding_dir /= sliding_dir.norm();

        Matrix<Real> slip_tangent_contribution(spatial_dimension,
                                               spatial_dimension);

        for (UInt alpha : arange(spatial_dimension - 1)) {
          Real n_T_alpha = sliding_dir.dot(_basis(alpha));
          Matrix<Real> e_outer_e(spatial_dimension, spatial_dimension);
          Matrix<Real> e_outer_n(e_outer_e);
          e_outer_n.outerProduct(_basis(alpha), _normal);
          // first we add the second term in Wriggers equation 8.45
          slip_tangent_contribution +=
              mu * this->penalty * n_T_alpha * e_outer_n;

          for (UInt beta : arange(spatial_dimension - 1)) {
            Real n_T_beta = sliding_dir.dot(_basis(beta));
            e_outer_e.outerProduct(_basis(alpha), _basis(beta));
            // then we add the first term
            slip_tangent_contribution +=
                this->friction_penalty *
                (Math::kronecker(alpha, beta) - n_T_alpha * n_T_beta) *
                e_outer_e;
          }
        }
        _tangent += slip_tangent_contribution;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */

INSTANTIATE_MATERIAL(cohesive_linear_friction_CTO,
                     MaterialCohesiveLinearFrictionCTO);

} // namespace akantu
