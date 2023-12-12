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
#include "dof_manager_default.hh"
#include "solver_vector_default.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_
#define AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline SparseSolverVectorArray::SparseSolverVectorArray(
    DOFManagerDefault & dof_manager, const ID & id)
    : SparseSolverVector(dof_manager, id) {}

/* -------------------------------------------------------------------------- */
inline SparseSolverVectorArray::SparseSolverVectorArray(
    const SparseSolverVectorArray & vector, const ID & id)
    : SparseSolverVector(vector, id) {}

/* -------------------------------------------------------------------------- */
template <class Array_>
SparseSolverVector &
SparseSolverVectorArrayTmpl<Array_>::operator+(const SparseSolverVector & y) {
  const auto & y_ = aka::as_type<SparseSolverVectorArray>(y);
  this->vector += y_.getVector();

  ++this->release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_>
SparseSolverVector &
SparseSolverVectorArrayTmpl<Array_>::copy(const SparseSolverVector & y) {
  const auto & y_ = aka::as_type<SparseSolverVectorArray>(y);
  this->vector.copy(y_.getVector());

  this->release_ = y.release();
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_>
inline Int SparseSolverVectorArrayTmpl<Array_>::size() const {
  return this->dof_manager.getSystemSize();
}

/* -------------------------------------------------------------------------- */
template <class Array_>
inline Int SparseSolverVectorArrayTmpl<Array_>::localSize() const {
  return dof_manager.getLocalSystemSize();
}

} // namespace akantu

#endif /* AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH_ */
