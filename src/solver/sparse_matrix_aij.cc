/**
 * @file   sparse_matrix_aij.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Aug 21 2015
 * @date last modification: Mon Dec 04 2017
 *
 * @brief  Implementation of the AIJ sparse matrix
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "sparse_matrix_aij.hh"
#include "aka_iterators.hh"
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
#include "solver_vector_default.hh"
#include "terms_to_assemble.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::SparseMatrixAIJ(const Communicator & communicator,
                                 const UInt m, const UInt n,
                                 const SizeType & size_type,
                                 const MatrixType & matrix_type, const ID & id)
    : SparseMatrix(communicator, m, n, size_type, matrix_type, id),
      irn(0, 1, id + ":irn"), jcn(0, 1, id + ":jcn"), a(0, 1, id + ":a") {
}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::SparseMatrixAIJ(const SparseMatrixAIJ & matrix, const ID & id)
    : SparseMatrix(matrix, id), irn(matrix.irn, id + ":irn"),
      jcn(matrix.jcn, id + ":jcn"), a(matrix.a, id + ":a") {}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::~SparseMatrixAIJ() = default;

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::saveProfile(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.open(filename.c_str());

  // write header
  if (communicator.whoAmI() == 0) {

    outfile << "%%MatrixMarket matrix coordinate pattern";
    if (this->matrix_type == _symmetric)
      outfile << " symmetric";
    else
      outfile << " general";
    outfile << std::endl;
    outfile << this->m << " " << this->m << " " << this->nb_non_zero
            << std::endl;
  }

  for (auto p : arange(communicator.getNbProc())) {
    // write content
    if (communicator.whoAmI() == p) {
      for (UInt i = 0; i < this->nb_non_zero; ++i) {
        outfile << this->irn.storage()[i] << " " << this->jcn.storage()[i]
                << " 1" << std::endl;
      }
    }
    communicator.barrier();
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();
  auto & comm = this->communicator;

  // open and set the properties of the stream
  std::ofstream outfile;

  if (0 == comm.whoAmI()) {
    outfile.open(filename.c_str());
  } else {
    outfile.open(filename.c_str(), std::ios_base::app);
  }

  outfile.precision(std::numeric_limits<Real>::digits10);
  // write header
  decltype(nb_non_zero) nnz = this->nb_non_zero;
  comm.allReduce(nnz);

  if (comm.whoAmI() == 0) {
    outfile << "%%MatrixMarket matrix coordinate real";
    if (this->matrix_type == _symmetric)
      outfile << " symmetric";
    else
      outfile << " general";
    outfile << std::endl;
    outfile << this->m << " " << this->n << " " << nnz << std::endl;
  }

  for (auto p : arange(comm.getNbProc())) {
    // write content
    if (comm.whoAmI() == p) {
      for (UInt i = 0; i < this->nb_non_zero; ++i) {
        outfile << this->irn(i) << " " << this->jcn(i) << " " << this->a(i)
                << std::endl;
      }
    }
    comm.barrier();
  }
  // time to end
  outfile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::copyContent(const SparseMatrix & matrix) {
  AKANTU_DEBUG_IN();
  const auto & mat = aka::as_type<SparseMatrixAIJ>(matrix);
  AKANTU_DEBUG_ASSERT(nb_non_zero == mat.getNbNonZero(),
                      "The to matrix don't have the same profiles");
  a = mat.getA();

  this->value_release++;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::copyProfile(const SparseMatrix & other) {
  auto & A = aka::as_type<SparseMatrixAIJ>(other);

  SparseMatrix::clearProfile();

  this->irn.copy(A.irn);
  this->jcn.copy(A.jcn);

  this->irn_jcn_k.clear();

  UInt i, j, k;
  for (auto && data : enumerate(irn, jcn)) {
    std::tie(k, i, j) = data;

    this->irn_jcn_k[this->key(i - 1, j - 1)] = k;
  }

  this->nb_non_zero = this->irn.size();
  this->a.resize(this->nb_non_zero);

  this->a.set(0.);
  this->m = A.m;
  this->n = A.n;

  this->profile_release = A.profile_release;
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
template <class MatrixType>
void SparseMatrixAIJ::addMeToTemplated(MatrixType & B, Real alpha) const {
  UInt i, j;
  Real A_ij;
  for (auto && tuple : zip(irn, jcn, a)) {
    std::tie(i, j, A_ij) = tuple;
    B.add(i - 1, j - 1, alpha * A_ij);
  }
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::addMeTo(SparseMatrix & B, Real alpha) const {

  if (aka::is_of_type<SparseMatrixAIJ>(B)) {
    this->addMeToTemplated<SparseMatrixAIJ>(aka::as_type<SparseMatrixAIJ>(B),
                                            alpha);
  } else {
    //    this->addMeToTemplated<SparseMatrix>(*this, alpha);
  }
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::mul(Real alpha) {
  this->a *= alpha;
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::clear() {
  a.set(0.);

  this->value_release++;
}

} // namespace akantu
