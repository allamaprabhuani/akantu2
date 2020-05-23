/**
 * @file   vector_petsc.cc
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
#include "vector_petsc.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
VectorPETSc::VectorPETSc(const Communicator & communicator, const UInt n,
                         const SizeType & size_type, const ID & id) {
  const auto & mpi_data =
      aka::as_type<MPICommunicatorData>(communicator.getCommunicatorData());
  mpi_comm = mpi_data.getMPICommunicator();

  PETSc_call(VecCreate, mpi_comm, &x);
  detail::PETScSetName(x, id);

  PETSc_call(VecSetFromOptions, x);

  switch (size_type) {
  case SizeType::_local: {
    PETSc_call(VecSetSizes, x, n, PETSC_DECIDE);
    n_local = n;
    PETSc_call(VecGetSize, x, &(this->n));
    break;
  }
  case SizeType::_global: {
    PETSc_call(VecSetSizes, x, PETSC_DECIDE, n);
    PETSc_call(VecGetLocalSize, x, &n_local);
    break;
  }
  }

  VecType vec_type;
  PETSc_call(VecGetType, x, &vec_type);
  if (std::string(vec_type) == std::string(VECMPI)) {
    // PetscInt lowest_gidx, highest_gidx;
    // PETSc_call(VecGetOwnershipRange, x, &lowest_gidx, &highest_gidx);

    // std::vector<PetscInt> ghost_idx;
    // for (auto && d : arange(local_system_size)) {
    //   int gidx = dof_manager.localToGlobalEquationNumber(d);
    //   if (gidx != -1) {
    //     if ((gidx < lowest_gidx) or (gidx >= highest_gidx)) {
    //       ghost_idx.push_back(gidx);
    //     }
    //   }
    // }

    // PETSc_call(VecMPISetGhost, x, ghost_idx.size(), ghost_idx.data());
  } else {
    std::vector<int> idx(n_local);
    std::iota(idx.begin(), idx.end(), 0);
    ISLocalToGlobalMapping is;
    PETSc_call(ISLocalToGlobalMappingCreate, PETSC_COMM_SELF, 1, idx.size(),
               idx.data(), PETSC_COPY_VALUES, &is);
    PETSc_call(VecSetLocalToGlobalMapping, x, is);
    PETSc_call(ISLocalToGlobalMappingDestroy, &is);
  }
}

/* -------------------------------------------------------------------------- */
VectorPETSc::VectorPETSc(const VectorPETSc & vector, const ID & id)
    : mpi_comm(vector.mpi_comm) {
  if (vector.x) {
    PETSc_call(VecDuplicate, vector.x, &x);
    PETSc_call(VecCopy, vector.x, x);
    detail::PETScSetName(x, id);
    n = vector.n;
    n_local = vector.n_local;
  }
}

/* -------------------------------------------------------------------------- */
void internal::VectorPETSc::print(const Vec & y, std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
  stream << space << "VectorPETSc [" << std::endl;
  //stream << space << " + id: " << id << std::endl;
  PETSc_call(PetscViewerPushFormat, PETSC_VIEWER_STDOUT_WORLD,
             PETSC_VIEWER_ASCII_INDEX);
  PETSc_call(VecView, y, PETSC_VIEWER_STDOUT_WORLD);
  PETSc_call(PetscViewerPopFormat, PETSC_VIEWER_STDOUT_WORLD);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
VectorPETSc::~VectorPETSc() {
  if (x) {
    PETSc_call(VecDestroy, &x);
  }
}

/* -------------------------------------------------------------------------- */
void VectorPETSc::clear() {
  PETSc_call(VecSet, x, 0.);
  applyModifications();
}

/* -------------------------------------------------------------------------- */
void VectorPETSc::applyModifications() {
  PETSc_call(VecAssemblyBegin, x);
  PETSc_call(VecAssemblyEnd, x);
  updateGhost();
}

/* -------------------------------------------------------------------------- */
void VectorPETSc::updateGhost() {
  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);
  if (x_ghosted) {
    PETSc_call(VecGhostUpdateBegin, x, INSERT_VALUES, SCATTER_FORWARD);
    PETSc_call(VecGhostUpdateEnd, x, INSERT_VALUES, SCATTER_FORWARD);
  }
  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
VectorPETSc::operator const Array<Real> &() const {
  const_cast<Array<Real> &>(this->cache).resize(local_size());

  auto xl = internal::make_petsc_local_vector(x);
  auto cachep = internal::make_petsc_wrapped_vector(this->cache);

  PETSc_call(VecCopy, cachep, xl);
  return cache;
}

/* -------------------------------------------------------------------------- */
VectorPETSc & VectorPETSc::operator=(const VectorPETSc & y) {
  if (size() != y.size()) {
    PETSc_call(VecDuplicate, y, &x);
  }

  PETSc_call(VecCopy, y.x, x);
  release_ = y.release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
// VectorPETSc & VectorPETSc::operator=(const VectorPETSc & y) {
//   const auto & y_ = aka::as_type<VectorPETSc>(y);
//   return operator=(y_);
// }

/* -------------------------------------------------------------------------- */
VectorPETSc & VectorPETSc::operator+(const VectorPETSc & y) {
  PETSc_call(VecAXPY, x, 1., y.x);
  release_ = y.release_;
  return *this;
}

} // namespace akantu
