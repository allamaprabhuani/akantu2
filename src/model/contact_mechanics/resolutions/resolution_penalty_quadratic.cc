/**
 * @file   resolution_penalty_quadratic.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Wed Jun 09 2021
 *
 * @brief  Specialization of the resolution class for the quadratic penalty
 * method
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
#include "resolution_penalty_quadratic.hh"
#include "element_class_helper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ResolutionPenaltyQuadratic::ResolutionPenaltyQuadratic(
    ContactMechanicsModel & model, const ID & id)
    : Parent(model, id) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::initialize() {}

/* -------------------------------------------------------------------------- */
Real ResolutionPenaltyQuadratic::computeNormalTraction(Real & gap) {
  return epsilon_n * (macaulay(gap) * macaulay(gap) + macaulay(gap));
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeNormalForce(
    const ContactElement & element, Vector<Real> & force) {
  force.zero();

  auto & gaps = model.getGaps();
  auto & projections = model.getProjections();
  auto & normals = model.getNormals();

  auto surface_dimension = spatial_dimension - 1;

  Real gap(gaps.begin()[element.slave]);
  Vector<Real> normal(normals.begin(spatial_dimension)[element.slave]);
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);

  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());

  // compute normal traction
  Real p_n = computeNormalTraction(gap);
  p_n *= nodal_area[element.slave];

  UInt nb_nodes_per_contact = element.getNbNodes();
  Matrix<Real> shape_matric(spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
  ResolutionUtils::computeShapeFunctionMatric(element, projection,
                                              shape_matric);

  force.mul<true>(shape_matric, normal, p_n);
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeTangentialForce(
    const ContactElement & element, Vector<Real> & force) {

  if (mu == 0) {
    return;
  }

  force.zero();

  UInt surface_dimension = spatial_dimension - 1;

  // compute tangents
  auto & projections = model.getProjections();
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);

  auto & normals = model.getNormals();
  Vector<Real> normal(normals.begin(spatial_dimension)[element.slave]);

  auto & tangents = model.getTangents();
  Matrix<Real> covariant_basis(
      tangents.begin(surface_dimension, spatial_dimension)[element.slave]);

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
  Vector<Real> tangential_traction(
      tangential_tractions.begin(surface_dimension)[element.slave]);
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
      Vector<Real> tmp(force.size());
      tmp.mul<true>(shape_matric, tangent_alpha, traction_beta);
      tmp *=
          contravariant_metric_tensor(alpha, beta) * nodal_area[element.slave];
      force += tmp;
    }
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeTangentialTraction(
    const ContactElement & element, const Matrix<Real> & covariant_basis,
    Vector<Real> & traction_tangential) {

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
  bool stick = traction_trial_norm <= mu * p_n;

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
void ResolutionPenaltyQuadratic::computeTrialTangentialTraction(
    const ContactElement & element, const Matrix<Real> & covariant_basis,
    Vector<Real> & traction) {

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

  auto current_tangent = covariant_basis.transpose();
  auto previous_tangent = previous_covariant_basis.transpose();

  for (auto alpha : arange(surface_dimension)) {
    Vector<Real> tangent_alpha(current_tangent(alpha));
    for (auto gamma : arange(surface_dimension)) {
      for (auto beta : arange(surface_dimension)) {
        Vector<Real> tangent_beta(previous_tangent(beta));
        auto t_alpha_t_beta = tangent_beta.dot(tangent_alpha);
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

  Vector<Real> temp(surface_dimension);
  temp.mul<false>(covariant_metric_tensor, increment_xi, epsilon_t);

  traction -= temp;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeStickTangentialTraction(
    const ContactElement & /*element*/, Vector<Real> & traction_trial,
    Vector<Real> & traction_tangential) {
  traction_tangential = traction_trial;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeSlipTangentialTraction(
    const ContactElement & element, const Matrix<Real> & covariant_basis,
    Vector<Real> & traction_trial, Vector<Real> & traction_tangential) {
  UInt surface_dimension = spatial_dimension - 1;

  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];

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

  auto slip_direction = traction_trial;
  slip_direction /= traction_trial_norm;

  Real p_n = computeNormalTraction(gap);
  traction_tangential = slip_direction;
  traction_tangential *= mu * p_n;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeNormalModuli(
    const ContactElement & element, Matrix<Real> & stiffness) {

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
  const ElementType & type = element.master.type;
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
  Matrix<Real> mat_n(normal.storage(), normal.size(), 1.);
  n_outer_n.mul<false, true>(mat_n, mat_n);

  Matrix<Real> tmp(spatial_dimension, spatial_dimension * nb_nodes_per_contact);
  tmp.mul<false, false>(n_outer_n, A);

  k_main.mul<true, false>(A, tmp);
  k_main *= epsilon_n * heaviside(gap) * (2 * gap + 1) * nodal_area;

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

    Matrix<Real> n_outer_t(spatial_dimension, spatial_dimension);
    Matrix<Real> mat_t(tangent.storage(), tangent.size(), 1.);
    n_outer_t.mul<false, true>(mat_n, mat_t);

    Matrix<Real> t_outer_n(spatial_dimension, spatial_dimension);
    t_outer_n.mul<false, true>(mat_t, mat_n);

    for (auto && values2 : enumerate(shape_derivatives.transpose())) {
      auto & beta = std::get<0>(values2);
      auto & dnds = std::get<1>(values2);
      // construct Aj from shape function wrt to jth natural
      // coordinate
      construct_Aj(dnds);

      Matrix<Real> tmp(spatial_dimension,
                       spatial_dimension * nb_nodes_per_contact);
      Matrix<Real> tmp1(nb_nodes_per_contact * spatial_dimension,
                        spatial_dimension * nb_nodes_per_contact);
      tmp.mul<false, false>(n_outer_t, A);
      tmp1.mul<true, false>(Aj, tmp);
      tmp1 *= contravariant_metric_tensor(alpha, beta);
      k_rot1 += tmp1;

      tmp.mul<false, false>(t_outer_n, Aj);
      tmp1.mul<true, false>(A, tmp);
      tmp1 *= contravariant_metric_tensor(alpha, beta);
      k_rot2 += tmp1;
    }
  }

  k_rot1 *= -epsilon_n * heaviside(gap) * (gap * gap + gap) * nodal_area;
  k_rot2 *= -epsilon_n * heaviside(gap) * (gap * gap + gap) * nodal_area;

  stiffness += k_main + k_rot1 + k_rot2;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeTangentialModuli(
    const ContactElement & element, Matrix<Real> & stiffness) {
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
void ResolutionPenaltyQuadratic::computeStickModuli(
    const ContactElement & element, Matrix<Real> & stiffness) {

  auto surface_dimension = spatial_dimension - 1;

  auto & projections = model.getProjections();
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);

  auto & nodal_areas = model.getNodalArea();
  auto & nodal_area = nodal_areas.begin()[element.slave];

  // method from Schweizerhof and A. Konyukhov, K. Schweizerhof
  // DOI 10.1007/s00466-004-0616-7 and DOI 10.1007/s00466-003-0515-3

  // construct A matrix
  const ElementType & type = element.master.type;
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

  for (auto && values1 : enumerate(covariant_basis.transpose())) {
    auto & alpha = std::get<0>(values1);
    auto & tangent_alpha = std::get<1>(values1);

    Matrix<Real> t_outer_t(spatial_dimension, spatial_dimension);
    Matrix<Real> mat_t_alpha(tangent_alpha.storage(), tangent_alpha.size(), 1.);

    for (auto && values2 : enumerate(covariant_basis.transpose())) {
      auto & beta = std::get<0>(values2);
      auto & tangent_beta = std::get<1>(values2);

      Matrix<Real> mat_t_beta(tangent_beta.storage(), tangent_beta.size(), 1.);
      t_outer_t.mul<false, true>(mat_t_alpha, mat_t_beta);

      Matrix<Real> tmp(spatial_dimension,
                       spatial_dimension * nb_nodes_per_contact);
      Matrix<Real> tmp1(nb_nodes_per_contact * spatial_dimension,
                        spatial_dimension * nb_nodes_per_contact);
      tmp.mul<false, false>(t_outer_t, A);
      tmp1.mul<true, false>(A, tmp);
      tmp1 *= contravariant_metric_tensor(alpha, beta);
      k_main += tmp1;
    }
  }

  k_main *= -epsilon_t;

  // construct 2nd part of the stick modulii
  auto & tangential_tractions = model.getTangentialTractions();
  Vector<Real> tangential_traction(
      tangential_tractions.begin(surface_dimension)[element.slave]);

  Matrix<Real> k_second(nb_nodes_per_contact * spatial_dimension,
                        nb_nodes_per_contact * spatial_dimension);

  for (auto alpha : arange(surface_dimension)) {

    Matrix<Real> k_sum(nb_nodes_per_contact * spatial_dimension,
                       nb_nodes_per_contact * spatial_dimension);

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
        Matrix<Real> mat_t_gamma(tangent_gamma.storage(), tangent_gamma.size(),
                                 1.);

        for (auto && values3 : enumerate(covariant_basis.transpose())) {
          auto & theta = std::get<0>(values3);
          auto & tangent_theta = std::get<1>(values3);

          Matrix<Real> mat_t_theta(tangent_theta.storage(),
                                   tangent_theta.size(), 1.);
          t_outer_t.mul<false, true>(mat_t_gamma, mat_t_theta);

          Matrix<Real> tmp(spatial_dimension,
                           spatial_dimension * nb_nodes_per_contact);
          Matrix<Real> tmp1(nb_nodes_per_contact * spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
          tmp.mul<false, false>(t_outer_t, Aj);
          tmp1.mul<true, false>(A, tmp);
          tmp1 *= contravariant_metric_tensor(alpha, theta) *
                  contravariant_metric_tensor(beta, gamma);

          Matrix<Real> tmp2(spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
          Matrix<Real> tmp3(nb_nodes_per_contact * spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
          tmp2.mul<false, false>(t_outer_t, A);
          tmp3.mul<true, false>(Aj, tmp2);
          tmp3 *= contravariant_metric_tensor(alpha, gamma) *
                  contravariant_metric_tensor(beta, theta);

          k_sum += tmp1 + tmp3;
        }
      }
    }

    k_second += tangential_traction[alpha] * k_sum;
  }

  stiffness += k_main * nodal_area - k_second * nodal_area;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::computeSlipModuli(
    const ContactElement & element, Matrix<Real> & stiffness) {

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

  // restructure normal as a matrix for an outer product
  Matrix<Real> mat_n(normal.storage(), normal.size(), 1.);

  // method from Schweizerhof and A. Konyukhov, K. Schweizerhof
  // DOI 10.1007/s00466-004-0616-7 and DOI 10.1007/s00466-003-0515-3

  // construct A matrix
  const ElementType & type = element.master.type;
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

  for (auto && values1 : enumerate(covariant_basis.transpose())) {
    auto & alpha = std::get<0>(values1);
    auto & tangent_alpha = std::get<1>(values1);

    Matrix<Real> mat_t_alpha(tangent_alpha.storage(), tangent_alpha.size(), 1.);

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
      Matrix<Real> mat_t_beta(tangent_beta.storage(), tangent_beta.size(), 1.);
      t_outer_n.mul<false, true>(mat_t_beta, mat_n);

      Matrix<Real> tmp(spatial_dimension,
                       spatial_dimension * nb_nodes_per_contact);
      Matrix<Real> tmp1(nb_nodes_per_contact * spatial_dimension,
                        spatial_dimension * nb_nodes_per_contact);
      tmp.mul<false, false>(t_outer_n, A);
      tmp1.mul<true, false>(A, tmp);

      tmp1 *= epsilon_n * mu * tangential_traction[alpha] *
              contravariant_metric_tensor(alpha, beta);
      tmp1 /= traction_norm;

      k_first += tmp1 * nodal_area;

      // eq 107a
      t_outer_t.mul<false, true>(mat_t_alpha, mat_t_beta);

      tmp.mul<false, false>(t_outer_t, A);
      tmp1.mul<true, false>(A, tmp);

      tmp1 *= epsilon_t * mu * p_n * contravariant_metric_tensor(alpha, beta);
      tmp1 /= traction_norm;

      k_second += tmp1 * nodal_area;

      for (auto && values3 : enumerate(covariant_basis.transpose())) {
        auto & gamma = std::get<0>(values3);
        auto & tangent_gamma = std::get<1>(values3);

        Matrix<Real> mat_t_gamma(tangent_gamma.storage(), tangent_gamma.size(),
                                 1.);

        for (auto && values4 : enumerate(covariant_basis.transpose())) {
          auto & theta = std::get<0>(values4);
          auto & tangent_theta = std::get<1>(values4);

          Matrix<Real> mat_t_theta(tangent_theta.storage(),
                                   tangent_theta.size(), 1.);
          t_outer_t.mul<false, true>(mat_t_gamma, mat_t_theta);

          // eq 107b
          tmp.mul<false, false>(t_outer_t, A);
          tmp1.mul<true, false>(A, tmp);

          tmp1 *= epsilon_t * mu * p_n * tangential_traction[alpha] *
                  tangential_traction[beta];
          tmp1 *= contravariant_metric_tensor(alpha, gamma) *
                  contravariant_metric_tensor(beta, theta);
          tmp1 /= pow(traction_norm, 3);

          k_third += tmp1 * nodal_area;

          // eq 107c
          tmp.mul<false, false>(t_outer_t, Aj);
          tmp1.mul<true, false>(A, tmp);
          tmp1 *= contravariant_metric_tensor(alpha, theta) *
                  contravariant_metric_tensor(beta, gamma);
          tmp1 *= mu * p_n * tangential_traction[alpha];
          tmp1 /= traction_norm;

          Matrix<Real> tmp2(spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
          Matrix<Real> tmp3(nb_nodes_per_contact * spatial_dimension,
                            spatial_dimension * nb_nodes_per_contact);
          tmp2.mul<false, false>(t_outer_t, A);
          tmp3.mul<true, false>(Aj, tmp2);
          tmp3 *= contravariant_metric_tensor(alpha, gamma) *
                  contravariant_metric_tensor(beta, theta);
          tmp3 *= mu * p_n * tangential_traction[alpha];
          tmp3 /= traction_norm;

          k_fourth += (tmp1 + tmp3) * nodal_area;
        }
      }
    }
  }

  stiffness += k_third + k_fourth - k_first - k_second;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::beforeSolveStep() {}

/* -------------------------------------------------------------------------- */
void ResolutionPenaltyQuadratic::afterSolveStep(
    __attribute__((unused)) bool converged) {}

INSTANTIATE_RESOLUTION(penalty_quadratic, ResolutionPenaltyQuadratic);

} // namespace akantu
