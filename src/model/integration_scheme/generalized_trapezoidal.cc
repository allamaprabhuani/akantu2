/**
 * @file   generalized_trapezoidal.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 04 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  implementation of inline functions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "generalized_trapezoidal.hh"
#include "mesh.hh"
#include "aka_array.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::predictor(Real delta_t, Array<Real> & u,
                                       Array<Real> & u_dot,
                                       const Array<bool> & blocked_dofs) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real * u_val = u.storage();
  Real * u_dot_val = u_dot.storage();
  bool * blocked_dofs_val = blocked_dofs.storage();

  for (UInt d = 0; d < nb_degree_of_freedom; d++) {
    if (!(*blocked_dofs_val)) {
      *u_val += (1. - alpha) * delta_t * *u_dot_val;
    }
    u_val++;
    u_dot_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::corrector(const SolutionType & type, Real delta_t,
                                       Array<Real> & u, Array<Real> & u_dot,
                                       const Array<bool> & blocked_dofs,
                                       const Array<Real> & delta) const {
  AKANTU_DEBUG_IN();

  switch (type) {
  case _temperature:
    this->allCorrector<_temperature>(delta_t, u, u_dot, blocked_dofs, delta);
    break;
  case _temperature_rate:
    this->allCorrector<_temperature_rate>(delta_t, u, u_dot, blocked_dofs,
                                          delta);
    break;
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
Real
GeneralizedTrapezoidal::getTemperatureCoefficient(const SolutionType & type,
                                                  Real delta_t) const {
  switch (type) {
  case _temperature:
    return 1.;
  case _temperature_rate:
    return alpha * delta_t;
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }
}

/* -------------------------------------------------------------------------- */
Real
GeneralizedTrapezoidal::getTemperatureRateCoefficient(const SolutionType & type,
                                                      Real delta_t) const {
  switch (type) {
  case _temperature:
    return 1. / (alpha * delta_t);
  case _temperature_rate:
    return 1.;
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }
}

/* -------------------------------------------------------------------------- */
template <IntegrationScheme::SolutionType type>
void GeneralizedTrapezoidal::allCorrector(Real delta_t, Array<Real> & u,
                                          Array<Real> & u_dot,
                                          const Array<bool> & blocked_dofs,
                                          const Array<Real> & delta) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real e = getTemperatureCoefficient(type, delta_t);
  Real d = getTemperatureRateCoefficient(type, delta_t);

  Real * u_val = u.storage();
  Real * u_dot_val = u_dot.storage();
  Real * delta_val = delta.storage();
  bool * blocked_dofs_val = blocked_dofs.storage();

  for (UInt dof = 0; dof < nb_degree_of_freedom; dof++) {
    if (!(*blocked_dofs_val)) {
      *u_val += e * *delta_val;
      *u_dot_val += d * *delta_val;
    }
    u_val++;
    u_dot_val++;
    delta_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::assembleJacobian(const SolutionType & type,
                                              Real delta_t) {
  AKANTU_DEBUG_IN();

  SparseMatrix & J = this->dof_manager.getMatrix("J");

  const SparseMatrix & M = this->dof_manager.getMatrix("M");
  const SparseMatrix & K = this->dof_manager.getMatrix("K");

  Real c = this->getTemperatureRateCoefficient(type, delta_t);
  Real e = this->getTemperatureCoefficient(type, delta_t);

  J.add(M, e);
  J.add(K, c);

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
