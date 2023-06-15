/**
 * @file   solver_vector_default_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 01 2019
 * @date last modification: Sat May 23 2020
 *
 * @brief  Solver vector interface to Array
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
#include "dof_manager_default.hh"
#include "solver_vector_default.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_
#define AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline SolverVectorArray::SolverVectorArray(DOFManagerDefault & dof_manager,
                                            const ID & id)
    : SolverVector(dof_manager, id) {}

/* -------------------------------------------------------------------------- */
inline SolverVectorArray::SolverVectorArray(const SolverVectorArray & vector,
                                            const ID & id)
    : SolverVector(vector, id) {}

/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector &
SolverVectorArrayTmpl<Array_>::operator+=(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  this->vector += y_.getVector();

  ++this->release_;
  return *this;
}
/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector &
SolverVectorArrayTmpl<Array_>::operator-=(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  this->vector -= y_.getVector();

  ++this->release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector &
SolverVectorArrayTmpl<Array_>::operator=(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  this->vector.copy(y_.getVector());

  this->release_ = y.release();
  return *this;
}
/* -------------------------------------------------------------------------- */
template <class Array_>
void SolverVectorArrayTmpl<Array_>::copy(const SolverVector & y) {
  this->operator=(y);
}

/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector & SolverVectorArrayTmpl<Array_>::operator*=(const Real & alpha) {
  this->vector *= alpha;

  ++this->release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_>
void SolverVectorArrayTmpl<Array_>::add(const SolverVector & y,
                                        const Real & alpha) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  for (auto && data : zip(this->vector, y_.getVector())) {
    auto && x = std::get<0>(data);
    auto && y = std::get<1>(data);
    x += y * alpha;
  }
  ++this->release_;
}

/* -------------------------------------------------------------------------- */
template <class Array_>
Real SolverVectorArrayTmpl<Array_>::dot(const SolverVector & y) const {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  Real sum = 0.;
  for (auto [d, a, b] : enumerate(this->vector, y_.getVector())) {
    sum += this->dof_manager.isLocalOrMasterDOF(d) * a * b;
  }
  return sum;
}

/* -------------------------------------------------------------------------- */
template <class Array_> Real SolverVectorArrayTmpl<Array_>::norm() const {
  auto a = this->dot(*this);
  return std::sqrt(a);
}

/* -------------------------------------------------------------------------- */
template <class Array_>
Real SolverVectorArrayTmpl<Array_>::normFreeDOFs() const {
  const auto & blocked_dofs = this->dof_manager.getBlockedDOFs();

  Real sum = 0.;
  for (auto [d, c, a] : enumerate(blocked_dofs, this->vector)) {
    bool is_local_node = this->dof_manager.isLocalOrMasterDOF(d);
    bool is_blocked_dof = c;
    if ((!is_blocked_dof) && is_local_node) {
      sum += a * a;
    }
  }
  return std::sqrt(sum);
}
/* -------------------------------------------------------------------------- */
template <class Array_> inline Int SolverVectorArrayTmpl<Array_>::size() const {
  return this->dof_manager.getSystemSize();
}

/* -------------------------------------------------------------------------- */
template <class Array_>
inline Int SolverVectorArrayTmpl<Array_>::localSize() const {
  return dof_manager.getLocalSystemSize();
}

} // namespace akantu

#endif /* AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_ */
