/**
 * @file   solver_vector_distributed.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed May 22 2024
 *
 * @brief  Solver vector interface for distributed arrays
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
#include "solver_vector_default.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_VECTOR_DISTRIBUTED_HH_
#define AKANTU_SOLVER_VECTOR_DISTRIBUTED_HH_

namespace akantu {

class SolverVectorDistributed : public SolverVectorDefault {
public:
  SolverVectorDistributed(DOFManagerDefault & dof_manager,
                          const ID & id = "solver_vector_mumps");

  SolverVectorDistributed(const SolverVectorDefault & vector,
                          const ID & id = "solver_vector_mumps");

  Array<Real> & getGlobalVector() override;
  void setGlobalVector(const Array<Real> & solution) override;

  bool isFinite() const override;

  Real dot(const SolverVector & y) const override;

  Real normFreeDOFs() const override;

  virtual bool isDistributed() const override { return true; }

protected:
  // full vector in case it needs to be centralized on master
  std::unique_ptr<Array<Real>> global_vector;
};

} // namespace akantu

#endif /* AKANTU_SOLVER_VECTOR_DISTRIBUTED_HH_ */
