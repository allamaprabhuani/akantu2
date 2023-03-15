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
Real ResolutionPenaltyQuadratic::computeNormalTraction(const Real & gap) const {
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
  auto A = ResolutionUtils::computeShapeFunctionMatrix(element, projection);

  // construct the main part of normal matrix
  Matrix<Real> k_main(A.cols(), A.cols());

  Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
  n_outer_n = normal * normal.transpose();

  k_main = (A.transpose() * n_outer_n * A) * epsilon_n * heaviside(gap) *
           (2 * gap + 1) * nodal_area;

  // construct the rotational part of the normal matrix
  auto & tangents = model.getTangents();
  auto && covariant_basis =
      tangents.begin(spatial_dimension, surface_dimension)[element.slave];

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
    // auto t_outer_n = tangent * normal.transpose();

    for (auto && [beta, Aj] : enumerate(Ajs)) {
      // construct Aj from shape function wrt to jth natural
      // coordinate
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
