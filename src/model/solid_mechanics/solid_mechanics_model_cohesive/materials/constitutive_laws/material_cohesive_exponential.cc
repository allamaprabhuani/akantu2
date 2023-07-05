/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_exponential.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialCohesiveExponential<dim>::MaterialCohesiveExponential(
    SolidMechanicsModel & model, const ID & id)
    : MaterialCohesive(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("beta", beta, Real(0.), _pat_parsable, "Beta parameter");

  this->registerParam("exponential_penalty", exp_penalty, true, _pat_parsable,
                      "Is contact penalty following the exponential law?");

  this->registerParam(
      "contact_tangent", contact_tangent, Real(1.0), _pat_parsable,
      "Ratio of contact tangent over the initial exponential tangent");

  // this->initInternalArray(delta_max, 1, _ek_cohesive);

  use_previous_delta_max = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialCohesiveExponential<dim>::initMaterial() {

  AKANTU_DEBUG_IN();
  MaterialCohesive::initMaterial();

  if ((exp_penalty) && (contact_tangent != 1)) {

    contact_tangent = 1;
    AKANTU_DEBUG_WARNING("The parsed paramter <contact_tangent> is forced to "
                         "1.0 when the contact penalty follows the exponential "
                         "law");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveExponential<dim>::computeTraction(ElementType el_type,
                                                       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// compute scalars
  Real beta2 = beta * beta;

  /// loop on each quadrature point
  for (auto && data : zip(make_view<dim>(tractions(el_type, ghost_type)),
                          make_view<dim>(opening(el_type, ghost_type)),
                          make_view<dim>(normals(el_type, ghost_type)),
                          delta_max(el_type, ghost_type),
                          delta_max.previous(el_type, ghost_type))) {
    auto & traction = std::get<0>(data);
    auto & opening = std::get<1>(data);
    auto & normal = std::get<2>(data);
    auto & delta_max = std::get<3>(data);
    auto & delta_max_prev = std::get<4>(data);

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening.dot(normal);
    auto && normal_opening = normal * normal_opening_norm;
    auto && tangential_opening = opening - normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    /**
     * compute effective opening displacement
     * @f$ \delta = \sqrt{
     * \beta^2 \Delta_t^2 + \Delta_n^2 } @f$
     */
    Real delta =
        std::sqrt(tangential_opening_norm * tangential_opening_norm * beta2 +
                  normal_opening_norm * normal_opening_norm);

    if ((normal_opening_norm < 0) &&
        (std::abs(normal_opening_norm) > Math::getTolerance())) {

      auto && delta_s = opening - normal_opening;
      delta = tangential_opening_norm * beta;

      computeCoupledTraction(traction, normal, delta, delta_s, delta_max,
                             delta_max_prev);

      computeCompressiveTraction(traction, normal, normal_opening_norm,
                                 opening);

    } else
      computeCoupledTraction(traction, normal, delta, opening, delta_max,
                             delta_max_prev);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class D1, class D2, class D3>
void MaterialCohesiveExponential<dim>::computeCoupledTraction(
    Eigen::MatrixBase<D1> & tract, const Eigen::MatrixBase<D2> & normal,
    Real delta, const Eigen::MatrixBase<D3> & opening, Real & delta_max_new,
    Real delta_max) {
  /// full damage case
  if (std::abs(delta) < Math::getTolerance()) {
    /// set traction to zero
    tract.zero();
    return;
  }

  /// element not fully damaged
  /**
   * Compute traction loading @f$ \mathbf{T} =
   * e \sigma_c \frac{\delta}{\delta_c} e^{-\delta/ \delta_c}@f$
   */
  /**
   * Compute traction unloading @f$ \mathbf{T} =
   *  \frac{t_{max}}{\delta_{max}} \delta @f$
   */
  Real beta2 = beta * beta;
  Real normal_open_norm = opening.dot(normal);

  delta_max_new = std::max(delta_max, delta);
  auto && op_n_n = normal * (1 - beta2) * normal_open_norm;
  tract = (beta2 * opening + op_n_n) * std::exp(1.) * sigma_c *
          std::exp(-delta_max_new / delta_c) / delta_c;
}

/* ------------------------------------------------------------------------- */
template <Int dim>
template <class D1, class D2, class D3>
void MaterialCohesiveExponential<dim>::computeCompressiveTraction(
    Eigen::MatrixBase<D1> & tract, const Eigen::MatrixBase<D2> & normal,
    Real delta_n, const Eigen::MatrixBase<D3> & /*opening*/) {
  Vector<Real> temp_tract(normal);

  if (exp_penalty) {
    temp_tract *= delta_n * std::exp(1) * sigma_c *
                  std::exp(-delta_n / delta_c) / delta_c;
  } else {
    Real initial_tg =
        contact_tangent * std::exp(1.) * sigma_c * delta_n / delta_c;
    temp_tract *= initial_tg;
  }

  tract += temp_tract;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveExponential<dim>::computeTangentTraction(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real beta2 = beta * beta;

  /**
   * compute tangent matrix  @f$ \frac{\partial \mathbf{t}}
   * {\partial \delta} = \hat{\mathbf{t}} \otimes
   * \frac{\partial (t/\delta)}{\partial \delta}
   * \frac{\hat{\mathbf{t}}}{\delta}+ \frac{t}{\delta}  [ \beta^2 \mathbf{I} +
   * (1-\beta^2) (\mathbf{n} \otimes \mathbf{n})] @f$
   **/

  /**
   * In which @f$
   *  \frac{\partial(t/ \delta)}{\partial \delta} =
   * \left\{\begin{array} {l l}
   *  -e  \frac{\sigma_c}{\delta_c^2  }e^{-\delta  /  \delta_c} &  \quad  if
   *  \delta \geq \delta_{max} \\
   *  0 & \quad if \delta < \delta_{max}, \delta_n > 0
   *  \end{array}\right. @f$
   **/

  for (auto && data : zip(make_view<dim, dim>(tangent_matrix),
                          make_view<dim>(opening(el_type, ghost_type)),
                          make_view<dim>(normals(el_type, ghost_type)),
                          delta_max.previous(el_type, ghost_type))) {
    auto && tangent = std::get<0>(data);
    auto && opening = std::get<1>(data);
    auto && normal = std::get<2>(data);
    auto && delta_max_prev = std::get<3>(data);

    Real normal_opening_norm = opening.dot(normal);

    auto && normal_opening = normal * normal_opening_norm;
    auto && tangential_opening = opening - normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    auto delta =
        std::sqrt(tangential_opening_norm * tangential_opening_norm * beta2 +
                  normal_opening_norm * normal_opening_norm);

    if ((normal_opening_norm < 0) &&
        (std::abs(normal_opening_norm) > Math::getTolerance())) {

      auto && delta_s = opening - normal_opening;
      delta = tangential_opening_norm * beta;

      computeCoupledTangent(tangent, normal, delta, delta_s, delta_max_prev);

      computeCompressivePenalty(tangent, normal, normal_opening_norm);
    } else
      computeCoupledTangent(tangent, normal, delta, opening, delta_max_prev);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class D1, class D2, class D3>
void MaterialCohesiveExponential<dim>::computeCoupledTangent(
    Eigen::MatrixBase<D1> & tangent, const Eigen::MatrixBase<D2> & normal,
    Real delta, const Eigen::MatrixBase<D3> & opening,
    Real /* delta_max_new */) {
  AKANTU_DEBUG_IN();

  auto beta2 = beta * beta;
  auto J = Matrix<Real, dim, dim>::Identity() * beta2;

  if (std::abs(delta) < Math::getTolerance()) {
    delta = Math::getTolerance();
  }

  auto opening_normal = opening.dot(normal);

  auto && delta_e = normal * opening_normal * (1. - beta2) + beta2 * opening;

  auto exponent = std::exp(1. - delta / delta_c) * sigma_c / delta_c;

  auto && first_term = normal * normal.transpose() * (1. - beta2) + J;

  auto && second_term = delta_e * delta_e.transpose() / delta / delta_c;

  tangent = (first_term - second_term) * exponent;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class D1, class D2>
void MaterialCohesiveExponential<dim>::computeCompressivePenalty(
    Eigen::MatrixBase<D1> & tangent, const Eigen::MatrixBase<D2> & normal,
    Real delta_n) {

  if (not exp_penalty) {
    delta_n = 0.;
  }

  auto && n_outer_n = normal * normal.transpose();

  auto normal_tg = contact_tangent * std::exp(1.) * sigma_c *
                   std::exp(-delta_n / delta_c) * (1. - delta_n / delta_c) /
                   delta_c;

  tangent = tangent + n_outer_n * normal_tg;
}

/* -------------------------------------------------------------------------- */
template class MaterialCohesiveExponential<1>;
template class MaterialCohesiveExponential<2>;
template class MaterialCohesiveExponential<3>;
static bool material_is_alocated_cohesive_exponential =
    instantiateMaterial<MaterialCohesiveExponential>("cohesive_exponential");

} // namespace akantu
