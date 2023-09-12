/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "resolution_penalty.hh"
#include "element_class_helper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ResolutionPenalty::ResolutionPenalty(ContactMechanicsModel & model,
                                     const ID & id)
    : Resolution(model, id) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::initialize() {
  this->registerParam("epsilon_n", epsilon_n, Real(0.),
                      _pat_parsable | _pat_modifiable,
                      "Normal penalty parameter");
  this->registerParam("epsilon_t", epsilon_t, Real(0.),
                      _pat_parsable | _pat_modifiable,
                      "Tangential penalty parameter");
}

/* -------------------------------------------------------------------------- */
Real ResolutionPenalty::computeNormalTraction(const Real & gap) const {
  return epsilon_n * macaulay(gap);
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalForce(const ContactElement & element,
                                           Vector<Real> & force) {
  const auto & gaps = model.getGaps();
  const auto & projections = model.getProjections();
  const auto & normals = model.getNormals();
  const auto & nodal_area = model.getNodalArea();

  auto surface_dimension = spatial_dimension - 1;

  auto gap = gaps(element.slave);
  auto && normal = normals.begin(spatial_dimension)[element.slave];
  auto && projection = projections.begin(surface_dimension)[element.slave];

  // compute normal traction
  auto p_n = computeNormalTraction(gap) * nodal_area[element.slave];

  auto shape_matrix =
      ResolutionUtils::computeShapeFunctionMatrix(element, projection);

  force = p_n * shape_matrix.transpose() * normal;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentialForce(const ContactElement & element,
                                               Vector<Real> & force) {
  force.zero();
  if (mu == 0) {
    return;
  }

  auto surface_dimension = spatial_dimension - 1;

  // compute covariant basis
  auto && projection =
      model.getProjections().begin(surface_dimension)[element.slave];

  auto && covariant_basis = model.getTangents().begin(
      spatial_dimension, surface_dimension)[element.slave];

  // check for no-contact to contact condition
  // need a better way to check if new node added is not presnt in the
  // previous master elemets
  auto & previous_master_elements = model.getPreviousMasterElements();
  if (element.slave >= previous_master_elements.size()) {
    return;
  }

  auto & previous_element = previous_master_elements[element.slave];
  if (previous_element.type == _not_defined) {
    return;
  }

  // compute tangential traction using return map algorithm
  auto & tangential_tractions = model.getTangentialTractions();
  auto && tangential_traction =
      tangential_tractions.begin(surface_dimension)[element.slave];
  this->computeTangentialTraction(element, covariant_basis,
                                  tangential_traction);

  auto shape_matrix =
      ResolutionUtils::computeShapeFunctionMatrix(element, projection);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  const auto & nodal_area = model.getNodalArea();

  for (auto && [alpha, tangent_alpha] : enumerate(covariant_basis)) {
    for (auto && [beta, traction_beta] : enumerate(tangential_traction)) {
      force += (traction_beta * shape_matrix.transpose() * tangent_alpha) *
               contravariant_metric_tensor(alpha, beta) *
               nodal_area[element.slave];
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename D>
void ResolutionPenalty::computeTangentialTraction(
    const ContactElement & element, const Matrix<Real> & covariant_basis,
    Eigen::MatrixBase<D> & traction_tangential) {

  Int surface_dimension = spatial_dimension - 1;

  const auto & gap = model.getGaps()(element.slave);

  // Return map algorithm is employed
  // compute trial traction
  Vector<Real> traction_trial(surface_dimension);
  this->computeTrialTangentialTraction(element, covariant_basis,
                                       traction_trial);

  // compute norm of trial traction
  Real traction_trial_norm = 0;
  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);
  for (auto i : arange(surface_dimension)) {
    for (auto j : arange(surface_dimension)) {
      traction_trial_norm += traction_trial[i] * traction_trial[j] *
                             contravariant_metric_tensor(i, j);
    }
  }
  traction_trial_norm = std::sqrt(traction_trial_norm);

  // check stick or slip condition
  auto & contact_state = model.getContactState();
  auto & state = contact_state.begin()[element.slave];

  Real p_n = computeNormalTraction(gap);
  bool stick = (traction_trial_norm <= mu * p_n);

  if (stick) {
    state = ContactState::_stick;
    computeStickTangentialTraction(element, traction_trial,
                                   traction_tangential);
  } else {
    state = ContactState::_slip;
    computeSlipTangentialTraction(element, covariant_basis, traction_trial,
                                  traction_tangential);
  }
}

/* -------------------------------------------------------------------------- */
template <typename D>
void ResolutionPenalty::computeTrialTangentialTraction(
    const ContactElement & element, const Matrix<Real> & current_tangent,
    Eigen::MatrixBase<D> & traction) {

  UInt surface_dimension = spatial_dimension - 1;

  auto & projections = model.getProjections();
  Vector<Real> current_projection(
      projections.begin(surface_dimension)[element.slave]);

  auto & previous_projections = model.getPreviousProjections();
  Vector<Real> previous_projection(
      previous_projections.begin(surface_dimension)[element.slave]);

  // method from Laursen et. al.
  /*auto covariant_metric_tensor =
  GeometryUtils::covariantMetricTensor(covariant_basis); auto
  increment_projection = current_projection - previous_projection;

  traction.mul<false>(covariant_metric_tensor, increment_projection, epsilon_t);

  auto & previous_tangential_tractions = model.getPreviousTangentialTractions();
  Vector<Real>
  previous_traction(previous_tangential_tractions.begin(surface_dimension)[element.slave]);
  traction = previous_traction + traction;*/

  // method from Schweizerhof
  auto covariant_metric_tensor =
      GeometryUtils::covariantMetricTensor(current_tangent);

  auto & previous_tangential_tractions = model.getPreviousTangentialTractions();
  auto && previous_traction(
      previous_tangential_tractions.begin(surface_dimension)[element.slave]);

  auto & previous_tangents = model.getPreviousTangents();
  auto && previous_tangent = previous_tangents.begin(
      spatial_dimension, surface_dimension)[element.slave];
  auto previous_contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(previous_tangent);

  for (auto alpha : arange(surface_dimension)) {
    for (auto gamma : arange(surface_dimension)) {
      for (auto beta : arange(surface_dimension)) {
        auto t_alpha_t_beta =
            previous_tangent(beta).dot(current_tangent(alpha));
        traction[alpha] += previous_traction[gamma] *
                           previous_contravariant_metric_tensor(gamma, beta) *
                           t_alpha_t_beta;
      }
    }
  }

  auto & previous_master_elements = model.getPreviousMasterElements();
  auto & previous_element = previous_master_elements[element.slave];

  auto coords = model.getMesh().extractNodalValuesFromElement(
      model.getContactDetector().getPositions(), previous_element);

  auto previous_real_projection = GeometryUtils::realProjection(
      coords, previous_element, previous_projection);

  auto current_real_projection =
      GeometryUtils::realProjection(coords, element.master, current_projection);

  auto increment_real = current_real_projection - previous_real_projection;
  Vector<Real> increment_xi(surface_dimension);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(current_tangent);

  // increment in natural coordinate
  for (auto beta : arange(surface_dimension)) {
    for (auto gamma : arange(surface_dimension)) {
      auto temp = increment_real.dot(current_tangent(gamma));
      temp *= contravariant_metric_tensor(beta, gamma);
      increment_xi[beta] += temp;
    }
  }

  traction -= epsilon_t * covariant_metric_tensor * increment_xi;
}

/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
void ResolutionPenalty::computeStickTangentialTraction(
    const ContactElement & /*element*/, Eigen::MatrixBase<D1> & traction_trial,
    Eigen::MatrixBase<D2> & traction_tangential) {
  traction_tangential = traction_trial;
}

/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
void ResolutionPenalty::computeSlipTangentialTraction(
    const ContactElement & element, const Matrix<Real> & covariant_basis,
    Eigen::MatrixBase<D1> & traction_trial,
    Eigen::MatrixBase<D2> & traction_tangential) {
  UInt surface_dimension = spatial_dimension - 1;

  auto & gap = model.getGaps()(element.slave);

  // compute norm of trial traction
  Real traction_trial_norm = 0;
  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  for (auto alpha : arange(surface_dimension)) {
    for (auto beta : arange(surface_dimension)) {
      traction_trial_norm += traction_trial[alpha] * traction_trial[beta] *
                             contravariant_metric_tensor(alpha, beta);
    }
  }
  traction_trial_norm = sqrt(traction_trial_norm);

  Real p_n = computeNormalTraction(gap);
  traction_tangential = traction_trial / traction_trial_norm * mu * p_n;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalModuli(const ContactElement & element,
                                            Matrix<Real> & stiffness) {

  auto surface_dimension = spatial_dimension - 1;

  auto & gaps = model.getGaps();
  Real gap(gaps.begin()[element.slave]);

  auto & projections = model.getProjections();
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);

  auto & nodal_areas = model.getNodalArea();
  auto & nodal_area = nodal_areas.begin()[element.slave];

  auto & normals = model.getNormals();
  Vector<Real> normal(normals.begin(spatial_dimension)[element.slave]);

  // method from Schweizerhof and A. Konyukhov, K. Schweizerhof
  // DOI 10.1007/s00466-004-0616-7 and DOI 10.1007/s00466-003-0515-3

  // construct A matrix
  auto A = ResolutionUtils::computeShapeFunctionMatrix(element, projection);

  // construct the main part of normal matrix
  Matrix<Real> k_main(A.cols(), A.cols());

  Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
  n_outer_n = normal * normal.transpose();

  k_main =
      (A.transpose() * n_outer_n * A) * epsilon_n * heaviside(gap) * nodal_area;

  // construct the rotational part of the normal matrix
  const auto & tangents = model.getTangents();
  auto && covariant_basis(
      tangents.begin(spatial_dimension, surface_dimension)[element.slave]);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  // consists of 2 rotational parts
  Matrix<Real> k_rot1(A.cols(), A.cols());
  Matrix<Real> k_rot2(A.cols(), A.cols());

  k_rot1.zero();
  k_rot2.zero();

  auto Ajs = ResolutionUtils::computeDerivativeShapeFunctionMatrix(element,
                                                                   projection);

  for (auto && [alpha, tangent] : enumerate(covariant_basis)) {
    auto n_outer_t = normal * tangent.transpose();
    auto t_outer_n = tangent * normal.transpose();

    for (auto && [beta, Aj] : enumerate(Ajs)) {
      // construct Aj from shape function wrt to jth natural
      // coordinate
      k_rot1 += (Aj.transpose() * n_outer_t * A) *
                contravariant_metric_tensor(alpha, beta);

      k_rot2 += (A.transpose() * t_outer_n * Aj) *
                contravariant_metric_tensor(alpha, beta);
    }
  }

  k_rot1 *= -epsilon_n * heaviside(gap) * gap * nodal_area;
  k_rot2 *= -epsilon_n * heaviside(gap) * gap * nodal_area;

  stiffness += k_main + k_rot1 + k_rot2;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentialModuli(const ContactElement & element,
                                                Matrix<Real> & stiffness) {
  if (mu == 0) {
    return;
  }

  stiffness.zero();

  auto & contact_state = model.getContactState();
  auto state = contact_state.begin()[element.slave];

  switch (state) {
  case ContactState::_stick: {
    computeStickModuli(element, stiffness);
    break;
  }
  case ContactState::_slip: {
    computeSlipModuli(element, stiffness);
    break;
  }
  default:
    break;
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeStickModuli(const ContactElement & element,
                                           Matrix<Real> & stiffness) {

  auto surface_dimension = spatial_dimension - 1;

  auto & projections = model.getProjections();
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);

  auto & nodal_areas = model.getNodalArea();
  auto & nodal_area = nodal_areas.begin()[element.slave];

  // method from Schweizerhof and A. Konyukhov, K. Schweizerhof
  // DOI 10.1007/s00466-004-0616-7 and DOI 10.1007/s00466-003-0515-3

  // construct A matrix
  auto A = ResolutionUtils::computeShapeFunctionMatrix(element, projection);

  // computing shape derivatives
  auto Ajs = ResolutionUtils::computeDerivativeShapeFunctionMatrix(element,
                                                                   projection);

  // tangents should have been calculated in normal modulii
  auto & tangents = model.getTangents();
  auto && covariant_basis(
      tangents.begin(spatial_dimension, surface_dimension)[element.slave]);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  // construct 1st part of the stick modulii
  Matrix<Real> k_main(A.cols(), A.cols());
  k_main.zero();

  for (auto && [alpha, tangent_alpha] : enumerate(covariant_basis)) {
    for (auto && [beta, tangent_beta] : enumerate(covariant_basis)) {
      auto t_outer_t = tangent_alpha * tangent_beta.transpose();

      k_main += (A.transpose() * t_outer_t * A) *
                contravariant_metric_tensor(alpha, beta);
    }
  }

  k_main *= -epsilon_t;

  // construct 2nd part of the stick modulii
  auto & tangential_tractions = model.getTangentialTractions();
  auto && tangential_traction =
      tangential_tractions.begin(surface_dimension)[element.slave];

  Matrix<Real> k_second(A.cols(), A.cols());
  k_second.zero();

  Matrix<Real> k_sum(A.cols(), A.cols());

  for (auto alpha : arange(surface_dimension)) {

    k_sum.zero();

    for (auto && [beta, Aj] : enumerate(Ajs)) {
      for (auto && [gamma, tangent_gamma] : enumerate(covariant_basis)) {
        Matrix<Real> t_outer_t(spatial_dimension, spatial_dimension);

        for (auto && values3 : enumerate(covariant_basis.transpose())) {
          auto & theta = std::get<0>(values3);
          auto & tangent_theta = std::get<1>(values3);

          t_outer_t = tangent_gamma * tangent_theta.transpose();

          k_sum += (A.transpose() * t_outer_t * Aj) *
                       contravariant_metric_tensor(alpha, theta) *
                       contravariant_metric_tensor(beta, gamma) +
                   (Aj.transpose() * t_outer_t * A) *
                       contravariant_metric_tensor(alpha, gamma) *
                       contravariant_metric_tensor(beta, theta);
        }
      }
    }

    k_second += tangential_traction[alpha] * k_sum;
  }

  stiffness += k_main * nodal_area - k_second * nodal_area;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeSlipModuli(const ContactElement & element,
                                          Matrix<Real> & stiffness) {

  auto surface_dimension = spatial_dimension - 1;

  auto & gaps = model.getGaps();
  Real gap(gaps.begin()[element.slave]);

  auto & nodal_areas = model.getNodalArea();
  auto & nodal_area = nodal_areas.begin()[element.slave];

  // compute normal traction
  Real p_n = computeNormalTraction(gap);

  auto & projections = model.getProjections();
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);

  auto & normals = model.getNormals();
  Vector<Real> normal(normals.begin(spatial_dimension)[element.slave]);

  // method from Schweizerhof and A. Konyukhov, K. Schweizerhof
  // DOI 10.1007/s00466-004-0616-7 and DOI 10.1007/s00466-003-0515-3

  // construct A matrix
  auto A = ResolutionUtils::computeShapeFunctionMatrix(element, projection);

  // computing shape derivatives
  auto Ajs = ResolutionUtils::computeDerivativeShapeFunctionMatrix(element,
                                                                   projection);

  // tangents should have been calculated in normal modulii
  auto && covariant_basis = model.getTangents().begin(
      spatial_dimension, surface_dimension)[element.slave];

  auto & tangential_tractions = model.getTangentialTractions();
  Vector<Real> tangential_traction(
      tangential_tractions.begin(surface_dimension)[element.slave]);

  // compute norm of trial traction
  Real traction_norm = 0;
  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  for (auto i : arange(surface_dimension)) {
    for (auto j : arange(surface_dimension)) {
      traction_norm += tangential_traction[i] * tangential_traction[j] *
                       contravariant_metric_tensor(i, j);
    }
  }
  traction_norm = sqrt(traction_norm);

  // construct four parts of stick modulii (eq 107,107a-c)
  Matrix<Real> k_first(A.cols(), A.cols());
  Matrix<Real> k_second(A.cols(), A.cols());
  Matrix<Real> k_third(A.cols(), A.cols());
  Matrix<Real> k_fourth(A.cols(), A.cols());

  k_first.zero();
  k_second.zero();
  k_first.zero();
  k_fourth.zero();

  Matrix<Real> t_outer_n(spatial_dimension, spatial_dimension);
  Matrix<Real> t_outer_t(spatial_dimension, spatial_dimension);

  for (auto && [alpha, tangent_alpha] : enumerate(covariant_basis)) {
    for (auto && [beta, tangent_beta, Aj] :
         zip(arange(surface_dimension), covariant_basis, Ajs)) {
      // eq 107
      t_outer_n = tangent_beta * normal.transpose();

      k_first += (A.transpose() * t_outer_n * A) * epsilon_n * mu *
                 tangential_traction[alpha] *
                 contravariant_metric_tensor(alpha, beta) / traction_norm *
                 nodal_area;

      // eq 107a
      t_outer_t = tangent_alpha * tangent_beta.transpose();

      k_second += (A.transpose() * t_outer_t * A) * epsilon_t * mu * p_n *
                  contravariant_metric_tensor(alpha, beta) / traction_norm *
                  nodal_area;

      for (auto && values3 : enumerate(covariant_basis)) {
        auto & gamma = std::get<0>(values3);
        auto & tangent_gamma = std::get<1>(values3);

        for (auto && values4 : enumerate(covariant_basis)) {
          auto & theta = std::get<0>(values4);
          auto & tangent_theta = std::get<1>(values4);

          t_outer_t = tangent_gamma * tangent_theta.transpose();

          // eq 107b
          k_third += (A.transpose() * t_outer_t * A) * epsilon_t * mu * p_n *
                     tangential_traction[alpha] * tangential_traction[beta] *
                     contravariant_metric_tensor(alpha, gamma) *
                     contravariant_metric_tensor(beta, theta) /
                     pow(traction_norm, 3) * nodal_area;

          // eq 107c
          k_fourth += ((A.transpose() * t_outer_t * Aj) *
                           contravariant_metric_tensor(alpha, theta) *
                           contravariant_metric_tensor(beta, gamma) +
                       (Aj.transpose() * t_outer_t * A) *
                           contravariant_metric_tensor(alpha, gamma) *
                           contravariant_metric_tensor(beta, theta)) *
                      nodal_area * mu * p_n * tangential_traction[alpha] /
                      traction_norm;
        }
      }
    }
  }

  stiffness += k_third + k_fourth - k_first - k_second;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::beforeSolveStep() {}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::afterSolveStep(__attribute__((unused)) bool converged) {
}

INSTANTIATE_RESOLUTION(penalty_linear, ResolutionPenalty);

} // namespace akantu
