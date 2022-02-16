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
Real ResolutionPenaltyQuadratic::computeNormalTraction(Real & gap) const {
  return epsilon_n * (macaulay(gap) * macaulay(gap) + macaulay(gap));
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

  k_main = (A.transpose() * n_outer_n * A) * epsilon_n * heaviside(gap) *
           (2 * gap + 1) * nodal_area;

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

  k_rot1 *= -epsilon_n * heaviside(gap) * (gap * gap + gap) * nodal_area;
  k_rot2 *= -epsilon_n * heaviside(gap) * (gap * gap + gap) * nodal_area;

  stiffness += k_main + k_rot1 + k_rot2;
}

INSTANTIATE_RESOLUTION(penalty_quadratic, ResolutionPenaltyQuadratic);

} // namespace akantu
