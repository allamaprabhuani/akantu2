/**
 * @file   resolution_penalty.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Wed Jun 09 2021
 *
 * @brief  Specialization of the resolution class for the penalty method
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
Real ResolutionPenalty::computeNormalTraction(Real & gap) const {
  return epsilon_n * macaulay(gap);
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalForce(const ContactElement & element,
                                           Vector<Real> & force) {

  force.zero();

  auto & gaps = model.getGaps();
  auto & projections = model.getProjections();
  auto & normals = model.getNormals();

  auto surface_dimension = spatial_dimension - 1;

  Real gap(gaps.begin()[element.slave]);
  auto && normal = normals.begin(spatial_dimension)[element.slave];
  auto && projection = projections.begin(surface_dimension)[element.slave];

  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());

  // compute normal traction
  Real p_n = computeNormalTraction(gap);
  p_n *= nodal_area[element.slave];

  UInt nb_nodes_per_contact = element.getNbNodes();
  Matrix<Real> shape_matric(spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
  ResolutionUtils::computeShapeFunctionMatric(element, projection,
                                              shape_matric);

  force = p_n * shape_matric.transpose() * normal;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentialForce(const ContactElement & element,
                                               Vector<Real> & force) {
  if (mu == 0) {
    return;
  }

  force.zero();

  auto surface_dimension = spatial_dimension - 1;

  // compute covariant basis
  auto & projections = model.getProjections();
  auto && projection = projections.begin(surface_dimension)[element.slave];

  // auto & normals = model.getNormals();
  // auto && normal = normals.begin(spatial_dimension)[element.slave];

  auto & tangents = model.getTangents();
  auto && covariant_basis =
      tangents.begin(surface_dimension, spatial_dimension)[element.slave];

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

  UInt nb_nodes_per_contact = element.getNbNodes();
  Matrix<Real> shape_matric(spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
  ResolutionUtils::computeShapeFunctionMatric(element, projection,
                                              shape_matric);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());

  for (auto && values1 : enumerate(covariant_basis.transpose())) {
    auto & alpha = std::get<0>(values1);
    auto & tangent_alpha = std::get<1>(values1);
    for (auto && values2 : enumerate(tangential_traction)) {
      auto & beta = std::get<0>(values2);
      auto & traction_beta = std::get<1>(values2);
      force += (traction_beta * shape_matric.transpose() * tangent_alpha) *
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

  UInt surface_dimension = spatial_dimension - 1;

  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];

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
  traction_trial_norm = sqrt(traction_trial_norm);

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
    const ContactElement & element, const Matrix<Real> & covariant_basis,
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
      GeometryUtils::covariantMetricTensor(covariant_basis);

  auto & previous_tangential_tractions = model.getPreviousTangentialTractions();
  Vector<Real> previous_traction(
      previous_tangential_tractions.begin(surface_dimension)[element.slave]);

  auto & previous_tangents = model.getPreviousTangents();
  Matrix<Real> previous_covariant_basis(previous_tangents.begin(
      surface_dimension, spatial_dimension)[element.slave]);
  auto previous_contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(previous_covariant_basis);

  auto current_tangent = covariant_basis.transpose().eval();
  auto previous_tangent = previous_covariant_basis.transpose().eval();

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

  Vector<Real> previous_real_projection(spatial_dimension);
  GeometryUtils::realProjection(
      model.getMesh(), model.getContactDetector().getPositions(),
      previous_element, previous_projection, previous_real_projection);

  Vector<Real> current_real_projection(spatial_dimension);
  GeometryUtils::realProjection(
      model.getMesh(), model.getContactDetector().getPositions(),
      element.master, current_projection, current_real_projection);

  auto increment_real = current_real_projection - previous_real_projection;
  Vector<Real> increment_xi(surface_dimension);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

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

  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];

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
  ElementType type = element.master.type;
  auto && shapes = ElementClassHelper<_ek_regular>::getN(projection, type);

  UInt nb_nodes_per_contact = element.getNbNodes();
  Matrix<Real> A(spatial_dimension, spatial_dimension * nb_nodes_per_contact);

  for (auto i : arange(nb_nodes_per_contact)) {
    for (auto j : arange(spatial_dimension)) {
      if (i == 0) {
        A(j, i * spatial_dimension + j) = 1;
        continue;
      }
      A(j, i * spatial_dimension + j) = -shapes[i - 1];
    }
  }

  // construct the main part of normal matrix
  Matrix<Real> k_main(nb_nodes_per_contact * spatial_dimension,
                      nb_nodes_per_contact * spatial_dimension);

  Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
  n_outer_n = normal * normal.transpose();

  k_main =
      (A.transpose() * n_outer_n * A) * epsilon_n * heaviside(gap) * nodal_area;

  // construct the rotational part of the normal matrix
  auto & tangents = model.getTangents();
  Matrix<Real> covariant_basis(
      tangents.begin(surface_dimension, spatial_dimension)[element.slave]);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  // computing shape derivatives
  auto && shape_derivatives =
      ElementClassHelper<_ek_regular>::getDNDS(projection, type);

  // consists of 2 rotational parts
  Matrix<Real> k_rot1(nb_nodes_per_contact * spatial_dimension,
                      nb_nodes_per_contact * spatial_dimension);
  Matrix<Real> k_rot2(nb_nodes_per_contact * spatial_dimension,
                      nb_nodes_per_contact * spatial_dimension);

  k_rot1.zero();
  k_rot2.zero();

  Matrix<Real> Aj(spatial_dimension, spatial_dimension * nb_nodes_per_contact);

  auto construct_Aj = [&](auto && dnds) {
    for (auto i : arange(nb_nodes_per_contact)) {
      for (auto j : arange(spatial_dimension)) {
        if (i == 0) {
          Aj(j, i * spatial_dimension + j) = 0;
          continue;
        }
        Aj(j, i * spatial_dimension + j) = dnds(i - 1);
      }
    }
  };

  for (auto && values1 : enumerate(covariant_basis.transpose())) {
    auto & alpha = std::get<0>(values1);
    auto & tangent = std::get<1>(values1);

    auto n_outer_t = normal * tangent.transpose();
    // auto t_outer_n = tangent * normal.transpose();

    for (auto && values2 : enumerate(shape_derivatives.transpose())) {
      auto & beta = std::get<0>(values2);
      auto & dnds = std::get<1>(values2);
      // construct Aj from shape function wrt to jth natural
      // coordinate
      construct_Aj(dnds);

      k_rot1 += (Aj.transpose() * n_outer_t * A) *
                contravariant_metric_tensor(alpha, beta);

      k_rot2 += (A.transpose() * n_outer_t * Aj) *
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
  ElementType type = element.master.type;
  auto && shapes = ElementClassHelper<_ek_regular>::getN(projection, type);

  UInt nb_nodes_per_contact = element.getNbNodes();
  Matrix<Real> A(spatial_dimension, spatial_dimension * nb_nodes_per_contact);

  for (auto i : arange(nb_nodes_per_contact)) {
    for (auto j : arange(spatial_dimension)) {
      if (i == 0) {
        A(j, i * spatial_dimension + j) = 1;
        continue;
      }
      A(j, i * spatial_dimension + j) = -shapes[i - 1];
    }
  }

  // computing shape derivatives
  auto && shape_derivatives =
      ElementClassHelper<_ek_regular>::getDNDS(projection, type);

  Matrix<Real> Aj(spatial_dimension, spatial_dimension * nb_nodes_per_contact);

  auto construct_Aj = [&](auto && dnds) {
    for (auto i : arange(nb_nodes_per_contact)) {
      for (auto j : arange(spatial_dimension)) {
        if (i == 0) {
          Aj(j, i * spatial_dimension + j) = 0;
          continue;
        }
        Aj(j, i * spatial_dimension + j) = dnds(i - 1);
      }
    }
  };

  // tangents should have been calculated in normal modulii
  auto & tangents = model.getTangents();
  Matrix<Real> covariant_basis(
      tangents.begin(surface_dimension, spatial_dimension)[element.slave]);

  auto contravariant_metric_tensor =
      GeometryUtils::contravariantMetricTensor(covariant_basis);

  // construct 1st part of the stick modulii
  Matrix<Real> k_main(nb_nodes_per_contact * spatial_dimension,
                      nb_nodes_per_contact * spatial_dimension);
  k_main.zero();

  for (auto && values1 : enumerate(covariant_basis.transpose())) {
    auto & alpha = std::get<0>(values1);
    auto & tangent_alpha = std::get<1>(values1);

    for (auto && values2 : enumerate(covariant_basis.transpose())) {
      auto & beta = std::get<0>(values2);
      auto & tangent_beta = std::get<1>(values2);

      auto t_outer_t = tangent_alpha * tangent_beta.transpose();

      k_main += (A.transpose() * t_outer_t * A) *
                contravariant_metric_tensor(alpha, beta);
    }
  }

  k_main *= -epsilon_t;

  // construct 2nd part of the stick modulii
  auto & tangential_tractions = model.getTangentialTractions();
  Vector<Real> tangential_traction(
      tangential_tractions.begin(surface_dimension)[element.slave]);

  Matrix<Real> k_second(nb_nodes_per_contact * spatial_dimension,
                        nb_nodes_per_contact * spatial_dimension);
  k_second.zero();

  for (auto alpha : arange(surface_dimension)) {

    Matrix<Real> k_sum(nb_nodes_per_contact * spatial_dimension,
                       nb_nodes_per_contact * spatial_dimension);
    k_sum.zero();

    for (auto && values1 : enumerate(shape_derivatives.transpose())) {
      auto & beta = std::get<0>(values1);
      auto & dnds = std::get<1>(values1);
      // construct Aj from shape function wrt to jth natural
      // coordinate
      construct_Aj(dnds);
      for (auto && values2 : enumerate(covariant_basis.transpose())) {
        auto & gamma = std::get<0>(values2);
        auto & tangent_gamma = std::get<1>(values2);

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
  ElementType type = element.master.type;
  auto && shapes = ElementClassHelper<_ek_regular>::getN(projection, type);

  UInt nb_nodes_per_contact = element.getNbNodes();
  Matrix<Real> A(spatial_dimension, spatial_dimension * nb_nodes_per_contact);

  for (auto i : arange(nb_nodes_per_contact)) {
    for (auto j : arange(spatial_dimension)) {
      if (i == 0) {
        A(j, i * spatial_dimension + j) = 1;
        continue;
      }
      A(j, i * spatial_dimension + j) = -shapes[i - 1];
    }
  }

  // computing shape derivatives
  auto && shape_derivatives =
      ElementClassHelper<_ek_regular>::getDNDS(projection, type);

  Matrix<Real> Aj(spatial_dimension, spatial_dimension * nb_nodes_per_contact);

  auto construct_Aj = [&](auto && dnds) {
    for (auto i : arange(nb_nodes_per_contact)) {
      for (auto j : arange(spatial_dimension)) {
        if (i == 0) {
          Aj(j, i * spatial_dimension + j) = 0;
          continue;
        }
        Aj(j, i * spatial_dimension + j) = dnds(i - 1);
      }
    }
  };

  // tangents should have been calculated in normal modulii
  auto & tangents = model.getTangents();
  Matrix<Real> covariant_basis(
      tangents.begin(surface_dimension, spatial_dimension)[element.slave]);

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
  Matrix<Real> k_first(nb_nodes_per_contact * spatial_dimension,
                       nb_nodes_per_contact * spatial_dimension);
  Matrix<Real> k_second(nb_nodes_per_contact * spatial_dimension,
                        nb_nodes_per_contact * spatial_dimension);
  Matrix<Real> k_third(nb_nodes_per_contact * spatial_dimension,
                       nb_nodes_per_contact * spatial_dimension);
  Matrix<Real> k_fourth(nb_nodes_per_contact * spatial_dimension,
                        nb_nodes_per_contact * spatial_dimension);

  k_first.zero();
  k_second.zero();
  k_first.zero();
  k_fourth.zero();

  for (auto && values1 : enumerate(covariant_basis.transpose())) {
    auto & alpha = std::get<0>(values1);
    auto & tangent_alpha = std::get<1>(values1);

    Matrix<Real> t_outer_n(spatial_dimension, spatial_dimension);
    Matrix<Real> t_outer_t(spatial_dimension, spatial_dimension);

    for (auto && values2 :
         zip(arange(surface_dimension), covariant_basis.transpose(),
             shape_derivatives.transpose())) {
      auto & beta = std::get<0>(values2);
      auto & tangent_beta = std::get<1>(values2);
      auto & dnds = std::get<2>(values2);
      // construct Aj from shape function wrt to jth natural
      // coordinate
      construct_Aj(dnds);

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

      for (auto && values3 : enumerate(covariant_basis.transpose())) {
        auto & gamma = std::get<0>(values3);
        auto & tangent_gamma = std::get<1>(values3);

        for (auto && values4 : enumerate(covariant_basis.transpose())) {
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
