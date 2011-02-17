/**
 * @file   sparse_matrix_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 21:09:38 2010
 *
 * @brief  implementation of inline methods of the SparseMatrix class
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
inline UInt SparseMatrix::addToProfile(UInt i, UInt j) {
  if(!irn_jcn_to_k) irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;

  std::pair<UInt, UInt> jcn_irn = std::make_pair(i, j);

  AKANTU_DEBUG_ASSERT(irn_jcn_to_k->find(jcn_irn) == irn_jcn_to_k->end(),
		      "Couple (i,j) = (" << i << "," << j << ") already in the profile");

  irn.push_back(i + 1);
  jcn.push_back(j + 1);

  //  a.resize(a.getSize() + 1);
  Real zero = 0;
  a.push_back(zero);

  (*irn_jcn_to_k)[jcn_irn] = nb_non_zero;

  nb_non_zero++;

  return nb_non_zero - 1;
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrix::clear() {
  memset(a.values, 0, nb_non_zero*sizeof(Real));
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrix::addToMatrix(UInt i, UInt j, Real value) {
  std::pair<UInt, UInt> jcn_irn = std::make_pair(i, j);
  std::map<std::pair<UInt, UInt>, UInt>::iterator irn_jcn_to_k_it = irn_jcn_to_k->find(jcn_irn);

  AKANTU_DEBUG_ASSERT(irn_jcn_to_k_it != irn_jcn_to_k->end(),
		      "Couple (i,j) = (" << i << "," << j << ") does not exist in the profile");

  a.values[irn_jcn_to_k_it->second] += value;
}

/* -------------------------------------------------------------------------- */
// inline void SparseMatrix::addToMatrixSym(const Vector<Real> & local_matrix,
// 					 const Element & element) {
//   AKANTU_DEBUG_ASSERT(element_to_sparse_profile[element.type] != NULL,
// 		      "No profile stored for this kind of element call first buildProfile()");

//   UInt nb_values_per_elem   = element_to_sparse_profile[element.type]->getNbComponent();
//   UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

//   Real * mat_val = local_matrix.values;
//   UInt * elem_to_sparse_val = element_to_sparse_profile[element.type]->values + element.element * nb_values_per_elem;
//   Real * a_val = a.values;

//   for (UInt j = 0; j < nb_nodes_per_element * nb_degre_of_freedom; ++j) {
//     UInt i_end = (sparse_matrix_type == _symmetric) ? j + 1 :
//       nb_nodes_per_element * nb_degre_of_freedom;
//     for (UInt i = 0; i < i_end; ++i) {
//       UInt k = *(elem_to_sparse_val++);
//       a_val[k] += *(mat_val++);
//     }
//   }
// }

/* -------------------------------------------------------------------------- */
inline void SparseMatrix::addToMatrix(Real * local_matrix,
				      const Element & element) {
  AKANTU_DEBUG_ASSERT(element_to_sparse_profile[element.type] != NULL,
		      "No profile stored for this kind of element call first buildProfile()");

  UInt nb_values_per_elem   = element_to_sparse_profile[element.type]->getNbComponent();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

  Real * mat_val = local_matrix;
  UInt * elem_to_sparse_val = element_to_sparse_profile[element.type]->values + element.element * nb_values_per_elem;
  Real * a_val = a.values;

  for (UInt j = 0; j < nb_nodes_per_element * nb_degre_of_freedom; ++j) {
    for (UInt i = 0; i < nb_nodes_per_element * nb_degre_of_freedom; ++i) {
      UInt k = *(elem_to_sparse_val++);
      a_val[k] += *(mat_val++);
    }
  }
}
