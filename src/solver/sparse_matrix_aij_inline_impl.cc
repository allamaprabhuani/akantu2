/**
 * @file   sparse_matrix_aij_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 18 00:42:45 2015
 *
 * @brief  Implementation of inline functions of SparseMatrixAIJ
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

__BEGIN_AKANTU__

inline UInt SparseMatrixAIJ::addToProfile(UInt i, UInt j) {
  KeyCOO jcn_irn = this->key(i, j);

  coordinate_list_map::iterator it = this->irn_jcn_k.find(jcn_irn);
  if (!(it == this->irn_jcn_k.end())) return it->second;

  this->irn.push_back(i + 1);
  this->jcn.push_back(j + 1);
  this->a.push_back(0.);

  this->irn_jcn_k[jcn_irn] = this->nb_non_zero;

  (this->nb_non_zero)++;

  return (this->nb_non_zero - 1);
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixAIJ::clearProfile() {
  SparseMatrix::clearProfile();

  this->irn_jcn_k.clear();

  this->irn.resize(0);
  this->jcn.resize(0);
  this->a.resize(0);
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixAIJ::addToMatrix(UInt i, UInt j, Real value) {
  KeyCOO jcn_irn = this->key(i, j);
  coordinate_list_map::iterator irn_jcn_k_it = this->irn_jcn_k.find(jcn_irn);

  AKANTU_DEBUG_ASSERT(irn_jcn_k_it != this->irn_jcn_k.end(),
                      "Couple (i,j) = (" << i << "," << j
                                         << ") does not exist in the profile");

  this->a(irn_jcn_k_it->second) += value;
}

/* -------------------------------------------------------------------------- */
inline Real SparseMatrixAIJ::operator()(UInt i, UInt j) const {
  KeyCOO jcn_irn = this->key(i, j);
  coordinate_list_map::const_iterator irn_jcn_k_it =
      this->irn_jcn_k.find(jcn_irn);

  if (irn_jcn_k_it == this->irn_jcn_k.end()) return 0.;
  return this->a(irn_jcn_k_it->second);
}

/* -------------------------------------------------------------------------- */
inline Real & SparseMatrixAIJ::operator()(UInt i, UInt j) {
  KeyCOO jcn_irn = this->key(i, j);
  coordinate_list_map::iterator irn_jcn_k_it = this->irn_jcn_k.find(jcn_irn);
  AKANTU_DEBUG_ASSERT(irn_jcn_k_it != this->irn_jcn_k.end(),
                      "Couple (i,j) = (" << i << "," << j
                                         << ") does not exist in the profile");

  return this->a(irn_jcn_k_it->second);
}

__END_AKANTU__
