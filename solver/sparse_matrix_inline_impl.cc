/**
 * @file   sparse_matrix_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 21:09:38 2010
 *
 * @brief  implementation of inline methods of the SparseMatrix class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
// UInt SparseMatrix::addToProfile(UInt i, UInt j) {
//   nb_non_zero++;
//   irn.push_back(i);
//   jcn.push_back(j);
//   a.resize(a.getSize() + 1); a.values[a.getSize() - 1] = 0;

//   return nb_non_zero;
// }

// /* -------------------------------------------------------------------------- */
// UInt SparseMatrix::addToProfileLocal(UInt i, UInt j) {
//   UInt gi = local_to_global->values[i];
//   UInt gj = local_to_global->values[j];
//   return addToProfile(gi, gj);
// }

/* -------------------------------------------------------------------------- */
inline void SparseMatrix::addToMatrix(const Vector<Real> & local_matrix,
 			       const Element & element) {
  UInt nb_values_per_elem   = element_to_sparse_profile[element.type]->getNbComponent();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

  Real * mat_val = local_matrix.values;
  UInt * elem_to_sparse_val = element_to_sparse_profile[element.type]->values + element.element * nb_values_per_elem;
  Real * a_val = a.values;

  for (UInt j = 0; j < nb_nodes_per_element * nb_degre_of_freedom; ++j) {
    UInt i_end = (sparse_matrix_type == _symmetric) ? j + 1 :
      nb_nodes_per_element * nb_degre_of_freedom;
    for (UInt i = 0; i < i_end; ++i) {
      UInt k = *(elem_to_sparse_val++);
      a_val[k] += *(mat_val++);
    }
  }
}
