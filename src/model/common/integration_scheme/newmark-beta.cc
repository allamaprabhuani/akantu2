/**
 * @file   newmark-beta.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 23 2015
 * @date last modification: Wed Mar 27 2019
 *
 * @brief  implementation of the  newmark-@f$\beta@f$ integration  scheme.  This
 * implementation is taken from Méthodes  numériques en mécanique des solides by
 * Alain Curnier \note{ISBN: 2-88074-247-1}
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
#include "newmark-beta.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NewmarkBeta::NewmarkBeta(DOFManager & dof_manager, const ID & dof_id,
                         Real alpha, Real beta)
    : IntegrationScheme2ndOrder(dof_manager, dof_id), beta(beta), alpha(alpha) {

  this->registerParam("alpha", this->alpha, alpha, _pat_parsmod,
                      "The alpha parameter");
  this->registerParam("beta", this->beta, beta, _pat_parsmod,
                      "The beta parameter");
}

/* -------------------------------------------------------------------------- */
/*
 * @f$ \tilde{u_{n+1}} = u_{n} +  \Delta t \dot{u}_n + \frac{\Delta t^2}{2}
 * \ddot{u}_n @f$
 * @f$ \tilde{\dot{u}_{n+1}} = \dot{u}_{n} +  \Delta t \ddot{u}_{n} @f$
 * @f$ \tilde{\ddot{u}_{n}} = \ddot{u}_{n} @f$
 */
void NewmarkBeta::predictor(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                            Array<Real> & u_dot_dot,
                            const Array<bool> & blocked_dofs) const {
  AKANTU_DEBUG_IN();

  auto nb_nodes = u.size();
  auto nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  auto * u_val = u.data();
  auto * u_dot_val = u_dot.data();
  auto * u_dot_dot_val = u_dot_dot.data();
  auto * blocked_dofs_val = blocked_dofs.data();

  for (Int d = 0; d < nb_degree_of_freedom; d++) {
    if (!(*blocked_dofs_val)) {
      Real dt_a_n = delta_t * *u_dot_dot_val;

      *u_val += (1 - k * alpha) * delta_t * *u_dot_val +
                (.5 - h * alpha * beta) * delta_t * dt_a_n;
      *u_dot_val = (1 - k) * *u_dot_val + (1 - h * beta) * dt_a_n;
      *u_dot_dot_val = (1 - h) * *u_dot_dot_val;
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NewmarkBeta::corrector(const SolutionType & type, Real delta_t,
                            Array<Real> & u, Array<Real> & u_dot,
                            Array<Real> & u_dot_dot,
                            const Array<bool> & blocked_dofs,
                            const Array<Real> & delta) const {
  AKANTU_DEBUG_IN();

  switch (type) {
  case _acceleration: {
    this->allCorrector<_acceleration>(delta_t, u, u_dot, u_dot_dot,
                                      blocked_dofs, delta);
    break;
  }
  case _velocity: {
    this->allCorrector<_velocity>(delta_t, u, u_dot, u_dot_dot, blocked_dofs,
                                  delta);
    break;
  }
  case _displacement: {
    this->allCorrector<_displacement>(delta_t, u, u_dot, u_dot_dot,
                                      blocked_dofs, delta);
    break;
  }
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real NewmarkBeta::getAccelerationCoefficient(const SolutionType & type,
                                             Real delta_t) const {
  switch (type) {
  case _acceleration:
    return 1.;
  case _velocity:
    return 1. / (beta * delta_t);
  case _displacement:
    return 1. / (alpha * beta * delta_t * delta_t);
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }
}

/* -------------------------------------------------------------------------- */
Real NewmarkBeta::getVelocityCoefficient(const SolutionType & type,
                                         Real delta_t) const {
  switch (type) {
  case _acceleration:
    return beta * delta_t;
  case _velocity:
    return 1.;
  case _displacement:
    return 1. / (alpha * delta_t);
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }
}

/* -------------------------------------------------------------------------- */
Real NewmarkBeta::getDisplacementCoefficient(const SolutionType & type,
                                             Real delta_t) const {
  switch (type) {
  case _acceleration:
    return alpha * beta * delta_t * delta_t;
  case _velocity:
    return alpha * delta_t;
  case _displacement:
    return 1.;
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }
}

/* -------------------------------------------------------------------------- */
template <IntegrationScheme::SolutionType type>
void NewmarkBeta::allCorrector(Real delta_t, Array<Real> & u,
                               Array<Real> & u_dot, Array<Real> & u_dot_dot,
                               const Array<bool> & blocked_dofs,
                               const Array<Real> & delta) const {
  AKANTU_DEBUG_IN();

  auto nb_nodes = u.size();
  auto nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  auto c = getAccelerationCoefficient(type, delta_t);
  auto d = getVelocityCoefficient(type, delta_t);
  auto e = getDisplacementCoefficient(type, delta_t);

  auto * u_val = u.data();
  auto * u_dot_val = u_dot.data();
  auto * u_dot_dot_val = u_dot_dot.data();
  auto * delta_val = delta.data();
  auto * blocked_dofs_val = blocked_dofs.data();

  for (Int dof = 0; dof < nb_degree_of_freedom; dof++) {
    if (!(*blocked_dofs_val)) {
      *u_val += e * *delta_val;
      *u_dot_val += d * *delta_val;
      *u_dot_dot_val += c * *delta_val;
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    delta_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NewmarkBeta::assembleJacobian(const SolutionType & type, Real delta_t) {
  AKANTU_DEBUG_IN();

  auto & J = this->dof_manager.getMatrix("J");

  const auto & M = this->dof_manager.getMatrix("M");

  auto c = this->getAccelerationCoefficient(type, delta_t);
  auto e = this->getDisplacementCoefficient(type, delta_t);

  bool does_j_need_update = false;
  does_j_need_update |= M.getRelease() != m_release;

  // in explicit this coefficient is exactly 0.
  if (not(e == 0.)) {
    const auto & K = this->dof_manager.getMatrix("K");
    does_j_need_update |= K.getRelease() != k_release;
  }

  if (this->dof_manager.hasMatrix("C")) {
    const auto & C = this->dof_manager.getMatrix("C");
    does_j_need_update |= C.getRelease() != c_release;
  }

  does_j_need_update |= this->dof_manager.hasBlockedDOFsChanged();

  if (not does_j_need_update) {
    AKANTU_DEBUG_OUT();
    return;
  }

  J.copyProfile(M);
  // J.zero();

  if (not(e == 0.)) {
    const auto & K = this->dof_manager.getMatrix("K");
    J.add(K, e);
    k_release = K.getRelease();
  }

  J.add(M, c);
  m_release = M.getRelease();

  if (this->dof_manager.hasMatrix("C")) {
    auto d = this->getVelocityCoefficient(type, delta_t);
    const auto & C = this->dof_manager.getMatrix("C");
    J.add(C, d);
    c_release = C.getRelease();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
