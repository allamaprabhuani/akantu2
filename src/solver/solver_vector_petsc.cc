/**
 * @file   solver_vector_petsc.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solver_vector_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
SolverVectorTmpl<VectorPETSc, DOFManagerPETSc>::SolverVectorTmpl(
    DOFManagerPETSc & dof_manager, const ID & id)
    : SolverVector(dof_manager, id),
      VectorPETSc(dof_manager.getCommunicator(),
                  dof_manager.getPureLocalSystemSize(), SizeType::_local, id),
      dof_manager(dof_manager) {
  auto local_system_size = dof_manager.getPureLocalSystemSize();
  VecType vec_type;
  PETSc_call(VecGetType, x, &vec_type);
  if (std::string(vec_type) == std::string(VECMPI)) {
    PetscInt lowest_gidx, highest_gidx;
    PETSc_call(VecGetOwnershipRange, x, &lowest_gidx, &highest_gidx);

    std::vector<PetscInt> ghost_idx;
    for (auto && d : arange(local_system_size)) {
      int gidx = dof_manager.localToGlobalEquationNumber(d);
      if (gidx != -1) {
        if ((gidx < lowest_gidx) or (gidx >= highest_gidx)) {
          ghost_idx.push_back(gidx);
        }
      }
    }

    PETSc_call(VecMPISetGhost, x, ghost_idx.size(), ghost_idx.data());
  }
}

/* -------------------------------------------------------------------------- */
template <>
void SolverVectorTmpl<VectorPETSc, DOFManagerPETSc>::resize() {
  // the arrays are destroyed and recreated in the dof manager
  // resize is so not implemented
  //AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
Int SolverVectorTmpl<VectorPETSc, DOFManagerPETSc>::localSize() const {
  return VectorPETSc::local_size();
}

/* -------------------------------------------------------------------------- */
template <>
void SolverVectorTmpl<VectorPETSc, DOFManagerPETSc>::setGlobalVector(const VectorPETSc &) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
VectorPETSc & SolverVectorTmpl<VectorPETSc, DOFManagerPETSc>::getGlobalVector() {
  AKANTU_TO_IMPLEMENT();
}

} // namespace akantu
