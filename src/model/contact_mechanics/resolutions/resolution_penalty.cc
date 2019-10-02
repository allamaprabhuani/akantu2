/**
 * @file   resolution_penalty.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Jan 7 2019
 * @date last modification: Mon Jan 7 2019
 *
 * @brief  Specialization of the resolution class for the penalty method
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "resolution_penalty.hh"

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
  this->registerParam("quadratic", quadratic, bool(false),
                      _pat_parsable | _pat_modifiable, "penalty function");
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalForce(Vector<Real> & force,
                                           Vector<Real> & n,
                                           ContactElement & element) {

  force.clear();
  Real tn = element.gap * epsilon_n;

  if (quadratic) {
    tn = macaulay(tn) * macaulay(element.gap);
  } else {
    tn = macaulay(tn);
  }
  for (UInt i : arange(force.size())) {
    force[i] += n[i] * tn;
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionalForce(Vector<Real> & force,
                                               Array<Real> & d_alpha,
                                               ContactElement & element) {

  Matrix<Real> m_alpha_beta(spatial_dimension - 1, spatial_dimension - 1);
  ResolutionUtils::computeMetricTensor(m_alpha_beta, element.tangents);
  computeFrictionalTraction(m_alpha_beta, element);

  auto & traction = element.traction;
  for (auto && values : zip(traction, make_view(d_alpha, d_alpha.size()))) {
    auto & t_s = std::get<0>(values);
    auto & d_s = std::get<1>(values);
    force += d_s * t_s;
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalModuli(Matrix<Real> & ke,
                                            Array<Real> & n_alpha,
                                            Array<Real> & d_alpha,
                                            Vector<Real> & n,
                                            ContactElement & element) {

  Real tn = element.gap * epsilon_n;
  tn = macaulay(tn);

  Matrix<Real> n_mat(n.storage(), n.size(), 1);
  ke.mul<false, true>(n_mat, n_mat);
  ke *= epsilon_n * heaviside(element.gap);
  
  /* for (auto && values : zip(make_view(n_alpha, n_alpha.size()),
                            make_view(d_alpha, d_alpha.size()))) {
    auto & n_s = std::get<0>(values);
    auto & d_s = std::get<1>(values);

    Matrix<Real> ns_mat(n_s.storage(), n_s.size(), 1);
    Matrix<Real> ds_mat(d_s.storage(), d_s.size(), 1);

    Matrix<Real> tmp1(n_s.size(), n_s.size());
    tmp1.mul<false, true>(ns_mat, ds_mat);

    Matrix<Real> tmp2(n_s.size(), n_s.size());
    tmp2.mul<false, true>(ds_mat, ns_mat);

    ke -= (tmp1 + tmp2) * tn;
  }

  auto surface_dimension = spatial_dimension - 1;
  Matrix<Real> m_alpha_beta(surface_dimension, surface_dimension);
  ResolutionUtils::computeMetricTensor(m_alpha_beta, element.tangents);

  for (auto && values :
       zip(arange(surface_dimension), make_view(n_alpha, n_alpha.size()))) {

    auto & s = std::get<0>(values);
    auto & n_s = std::get<1>(values);

    Matrix<Real> ns_mat(n_s.storage(), n_s.size(), 1);

    for (auto && tuple :
         zip(arange(surface_dimension), make_view(n_alpha, n_alpha.size()))) {
      auto & t = std::get<0>(tuple);
      auto & n_t = std::get<1>(tuple);

      Matrix<Real> nt_mat(n_t.storage(), n_t.size(), 1);

      Matrix<Real> tmp3(n_s.size(), n_s.size());
      tmp3.mul<false, true>(ns_mat, nt_mat);

      tmp3 *= m_alpha_beta(s, t) * tn * element.gap;
      ke += tmp3;
    }
    }*/
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionalModuli(
    Matrix<Real> & /*ke*/, Array<Real> & t_alpha_beta,
    Array<Real> & n_alpha_beta, Array<Real> & /*n_alpha*/,
    Array<Real> & d_alpha, Matrix<Real> & phi, Vector<Real> & n,
    ContactElement & element) {

  auto k_common =
      computeCommonModuli(t_alpha_beta, n_alpha_beta, d_alpha, n, element);

  const auto & type = element.master.type;
  const auto & conn = element.connectivity;

  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;

  Matrix<Real> m_alpha_beta(surface_dimension, surface_dimension);
  ResolutionUtils::computeMetricTensor(m_alpha_beta, element.tangents);

  Array<Real> g_alpha(conn.size() * spatial_dimension, surface_dimension);
  ResolutionUtils::computeGalpha(g_alpha, t_alpha_beta, d_alpha, phi, element);

  /*Matrix<Real> k_t;
  bool stick = computeFrictionalTraction(m_alpha_beta, element);

  if(stick)
    k_t = computeStickModuli(g_alpha, d_alpha, m_alpha_beta);
  else
    k_t = computeSlipModuli(g_alpha, d_alpha, m_alpha_beta, element);*/
}

/* -------------------------------------------------------------------------- */
bool ResolutionPenalty::computeFrictionalTraction(Matrix<Real> & m_alpha_beta,
                                                  ContactElement & element) {

  Real tn = element.gap * epsilon_n;
  tn = macaulay(tn);

  auto delta_xi = element.projection - element.previous_projection;

  Vector<Real> trial_traction(delta_xi.size());

  trial_traction.mul<false>(m_alpha_beta, delta_xi, epsilon_n);
  trial_traction += element.traction;

  auto trial_slip_function = trial_traction.norm() - mu * tn;

  bool stick = false;
  if (trial_slip_function <= 0) {
    element.traction = trial_traction;
    stick = true;
  } else {
    element.traction = mu * tn * trial_traction / trial_traction.norm();
  }

  return stick;
}

/* -------------------------------------------------------------------------- */
Array<Real> ResolutionPenalty::computeCommonModuli(
    Array<Real> & /*t_alpha_beta*/, Array<Real> & /*n_alpha_beta*/,
    Array<Real> & d_alpha, Vector<Real> & /*n*/, ContactElement & /*element*/) {

  Array<Real> kt_alpha(spatial_dimension - 1, d_alpha.size() * d_alpha.size(),
                       "k_T_alpha");

  // auto t_alpha_beta_size = t_alpha_beta.size() * (spatial_dimension - 1);

  // auto & tangents = element.tangents;
  /*for(auto && values :
    zip(tangents.transpose(),
        make_view(kt_alpha, kt_alpha.size()),
        make_view(t_alpha_beta, t_alpha_beta_size),
        make_view(n_alpha_beta, n_alpha_beta_size))) {
    auto & tangent_s = std::get<0>(values);
    auto & kt_s      = std::get<1>(values);
    auto & t_alpha_s = std::get<2>(values);
    auto & n_alpha_s = std::get<3>(values);

    Matrix<Real> kt_s_mat(kt_s.storage(), d_alpha.size(), d_alpha.size());

    // loop over beta
    for(auto && tuple :
      zip(make_view(d_alpha, d_alpha.size()),
          make_view(n_alpha_ ))) {
      auto & d_s = std::get<0>(tuple);
      Matrix<Real> tmp(d_s.size(), d_s.size());

      // loop over gamma
      for(auto && entry :
        make_view(d_alpha, d_alpha.size())) {
    auto & d_t = std::get<0>(entry);


    // compute constant
    Matrix<Real> tmp2(d_t.size(), d_t.size());
    tmp2.mul<false, true>(d_s, d_t);
    kt_s_mat += tmp2;
      }
      }
      }*/

  return kt_alpha;
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ResolutionPenalty::computeStickModuli(
    Array<Real> & g_alpha, Array<Real> & d_alpha, Matrix<Real> & m_alpha_beta) {

  Matrix<Real> k_stick(d_alpha.size(), d_alpha.size());

  for (auto && values :
       zip(arange(d_alpha.getNbComponent()), make_view(d_alpha, d_alpha.size()),
           make_view(g_alpha, g_alpha.size()))) {
    auto & s = std::get<0>(values);
    auto & d_s = std::get<1>(values);
    auto & g_s = std::get<2>(values);

    Matrix<Real> ds_mat(d_s.storage(), d_s.size(), 1);
    Matrix<Real> gs_mat(g_s.storage(), g_s.size(), 1);

    Matrix<Real> tmp1(d_s.size(), d_s.size());
    tmp1.mul<false, true>(ds_mat, gs_mat);

    k_stick += tmp1;

    for (auto && tuple : enumerate(make_view(d_alpha, d_alpha.size()))) {
      auto & t = std::get<0>(tuple);
      auto & d_t = std::get<1>(tuple);

      Matrix<Real> dt_mat(d_t.storage(), d_t.size(), 1);
      Matrix<Real> tmp2(d_t.size(), d_t.size());

      tmp2.mul<false, true>(ds_mat, dt_mat);
      k_stick += tmp2 * m_alpha_beta(s, t);
    }
  }

  k_stick *= epsilon_t;
  return k_stick;
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ResolutionPenalty::computeSlipModuli(
						  Array<Real> & /*g_alpha*/, Array<Real> & d_alpha,
						  Matrix<Real> & /*m_alpha_beta*/, ContactElement & /*element*/) {

  /*Real tn = element.gap * epsilon;
  tn = macaulay(tn);

  Real factor;
  factor = epsilon_t * mu * tn;

  auto p_t = element.traction;
  p_t /= p_t.norm();*/

  Matrix<Real> k_slip(d_alpha.size(), d_alpha.size());

  /*
  // loop for alpha
  for(auto && value :
    make_view(d_alpha, d_alpha.size())) {
    auto & d_s = std::get<0>(value);

    // loop for beta
    for(auto && tuple :
      zip(arange(spatial_dimension - 1),
          make_view(d_alpha, d_alpha.size()),
          make_view(g_alpha, g_alpha.size()))) {
      auto & beta = std::get<0>(tuple);
      auto & d_beta = std::get<1>(tuple);
      auto & g_beta = std::get<2>(tuple);

      // loop for gamma
      for(auto && entry :
        zip(arange(spatial_dimension - 1),
        make_view(d_alpha, d_alpha.size()))) {
    auto & gamma = std::get<0>(entry);
    auto & d_gamma = std::get<1>(entry);

      }
    }


    }*/

  return k_slip;
}

INSTANTIATE_RESOLUTION(penalty, ResolutionPenalty);

} // namespace akantu
