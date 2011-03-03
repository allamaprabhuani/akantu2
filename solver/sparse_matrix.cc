/**
 * @file   sparse_matrix.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Oct 26 18:25:07 2010
 *
 * @brief  implementation of the SparseMatrix class
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
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "sparse_matrix.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(const Mesh & mesh,
			   const SparseMatrixType & sparse_matrix_type,
			   UInt nb_degre_of_freedom,
			   const SparseMatrixID & id,
			   const MemoryID & memory_id) :
  Memory(memory_id), id(id),
  sparse_matrix_type(sparse_matrix_type),
  nb_degre_of_freedom(nb_degre_of_freedom),
  mesh(&mesh),
  nb_non_zero(0),
  irn(0,1,"irn"), jcn(0,1,"jcn"), a(0,1,"A"),
  irn_jcn_to_k(NULL) {
  AKANTU_DEBUG_IN();

  size = mesh.getNbGlobalNodes()*nb_degre_of_freedom;

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->element_to_sparse_profile[t] = NULL;
  }

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  nb_proc = comm->getNbProc();

  irn_save = NULL;
  jcn_save = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(UInt size,
			   const SparseMatrixType & sparse_matrix_type,
			   UInt nb_degre_of_freedom,
			   const SparseMatrixID & id,
			   const MemoryID & memory_id) :
  Memory(memory_id), id(id),
  sparse_matrix_type(sparse_matrix_type),
  nb_degre_of_freedom(nb_degre_of_freedom),
  mesh(NULL), size(size*nb_degre_of_freedom),
  nb_non_zero(0),
  irn(0,1,"irn"), jcn(0,1,"jcn"), a(0,1,"A"),
  irn_jcn_to_k(NULL) {
  AKANTU_DEBUG_IN();

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->element_to_sparse_profile[t] = NULL;
  }

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  nb_proc = comm->getNbProc();

  irn_save = NULL;
  jcn_save = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(const SparseMatrix & matrix) :
  Memory(matrix.getMemoryID()), sparse_matrix_type(matrix.getSparseMatrixType()),
  nb_degre_of_freedom(matrix.getNbDegreOfFreedom()),
  mesh(NULL), size(matrix.getSize()), nb_non_zero(matrix.getNbNonZero()),
  irn(matrix.getIRN(), true), jcn(matrix.getJCN(), true), a(matrix.getA(), true),
  irn_jcn_to_k(NULL) {
  AKANTU_DEBUG_IN();

  irn_save = NULL;
  jcn_save = NULL;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
SparseMatrix::~SparseMatrix() {
  AKANTU_DEBUG_IN();

  if (irn_jcn_to_k) delete irn_jcn_to_k;

  if(irn_save) delete irn_save;
  if(jcn_save) delete jcn_save;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::buildProfile() {
  AKANTU_DEBUG_IN();

  if(irn_jcn_to_k) delete irn_jcn_to_k;
  irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;

  //  std::map<std::pair<UInt, UInt>, UInt> irn_jcn_to_k;
  std::map<std::pair<UInt, UInt>, UInt>::iterator irn_jcn_to_k_it;

  nb_non_zero = 0;
  irn.resize(0);
  jcn.resize(0);

  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if (Mesh::getSpatialDimension(*it) != mesh->getSpatialDimension()) continue;

    UInt nb_element = mesh->getNbElement(*it);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

    if (element_to_sparse_profile[*it] == NULL) {
      std::stringstream sstr;

      UInt nb_values_per_elem = nb_degre_of_freedom * nb_nodes_per_element;
      // if (sparse_matrix_type == _symmetric) {
      // 	nb_values_per_elem = (nb_values_per_elem * (nb_values_per_elem + 1)) / 2;
      // } else {
      nb_values_per_elem *= nb_values_per_elem;
      // }
      sstr << id << ":" << "element_to_sparse_profile:" << *it;
      element_to_sparse_profile[*it] = &(alloc<UInt>(sstr.str(),
    						     nb_element,
    						     nb_values_per_elem));
    }

    UInt * global_nodes_ids_val = NULL;
    if(nb_proc > 1) global_nodes_ids_val = mesh->getGlobalNodesIds().values;
    UInt * elem_to_sparse_val = element_to_sparse_profile[*it]->values;
    UInt * conn_val = mesh->getConnectivity(*it).values;

    for (UInt e = 0; e < nb_element; ++e) {                                     // loop on element
      UInt count = 0;
      for (UInt j = 0; j < nb_nodes_per_element; ++j) {                         // loop on local column
	UInt n_j = (nb_proc == 1) ? conn_val[j] : global_nodes_ids_val[conn_val[j]];
	UInt c_jcn = n_j * nb_degre_of_freedom;

	for (UInt d_j = 0; d_j < nb_degre_of_freedom; ++d_j, ++c_jcn) {         // loop on degre of freedom
	  //	  UInt i_end = (sparse_matrix_type == _symmetric) ? j + 1 : nb_nodes_per_element;
	  UInt i_end = nb_nodes_per_element;

	  for (UInt i = 0; i < i_end; ++i) {                                    // loop on rows
	    UInt n_i = (nb_proc == 1) ? conn_val[i] : global_nodes_ids_val[conn_val[i]];
	    UInt c_irn = n_i * nb_degre_of_freedom;
	    //	    UInt d_i_end = (sparse_matrix_type == _symmetric && i == j) ? d_j + 1 : nb_degre_of_freedom;
	    UInt d_i_end = nb_degre_of_freedom;

	    for (UInt d_i = 0; d_i < d_i_end; ++d_i, ++c_irn) {                 // loop on degre of freedom

	      std::pair<UInt, UInt> jcn_irn;
	      if ((sparse_matrix_type == _symmetric) && c_irn < c_jcn) jcn_irn = std::make_pair(c_jcn, c_irn);
	      else jcn_irn = std::make_pair(c_irn, c_jcn);
	      irn_jcn_to_k_it = irn_jcn_to_k->find(jcn_irn);

	      if (irn_jcn_to_k_it == irn_jcn_to_k->end()) {
		*elem_to_sparse_val++ = nb_non_zero;
		count++;
		(*irn_jcn_to_k)[jcn_irn] = nb_non_zero;
		irn.push_back(c_irn + 1);
		jcn.push_back(c_jcn + 1);
		nb_non_zero++;
	      } else {
		*elem_to_sparse_val++ = irn_jcn_to_k_it->second;
		count++;
	      }
	    }
	  }
	}
      }
      conn_val += nb_nodes_per_element;
    }
  }

  a.resize(nb_non_zero);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::applyBoundary(const Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  Int * irn_val = irn.values;
  Int * jcn_val = jcn.values;
  Real * a_val   = a.values;

  for (UInt i = 0; i < nb_non_zero; ++i) {
    if(boundary.values[*irn_val - 1] || boundary.values[*jcn_val - 1]) {
      *a_val *= (*irn_val == *jcn_val);
    }
    irn_val++; jcn_val++; a_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::removeBoundary(const Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  if(irn_save) delete irn_save;
  if(jcn_save) delete jcn_save;

  irn_save = new Vector<Int>(irn, true);
  jcn_save = new Vector<Int>(jcn, true);

  UInt n = boundary.getSize()*boundary.getNbComponent();

  UInt * perm = new UInt[n];

  size_save = size;
  size = 0;
  for (UInt i = 0; i < n; ++i) {
    if(!boundary.values[i]) {
      perm[i] = size;
      std::cout <<  "perm["<< i <<"] = " << size << std::endl;
      size++;
    }
  }

  for (UInt i = 0; i < nb_non_zero;) {
    if(boundary.values[irn.at(i) - 1] || boundary.values[jcn.at(i) - 1]) {
      irn.erase(i);
      jcn.erase(i);
      a.erase(i);
      nb_non_zero--;
    } else {
      irn.values[i] = perm[irn.values[i] - 1] + 1;
      jcn.values[i] = perm[jcn.values[i] - 1] + 1;
      i++;
    }
  }

  delete [] perm;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::restoreProfile() {
  AKANTU_DEBUG_IN();

  irn.resize(irn_save->getSize());
  jcn.resize(jcn_save->getSize());

  nb_non_zero = irn.getSize();
  a.resize(nb_non_zero);
  size = size_save;

  memcpy(irn.values, irn_save->values, irn.getSize()*sizeof(Int));
  memcpy(jcn.values, jcn_save->values, jcn.getSize()*sizeof(Int));

  delete irn_save; irn_save = NULL;
  delete jcn_save; jcn_save = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::saveProfile(const std::string & filename) {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.open(filename.c_str());

  outfile << "%%MatrixMarket matrix coordinate pattern";

  if(sparse_matrix_type == _symmetric) outfile << " symmetric";
  else outfile << " general";
  outfile << std::endl;

  UInt m = size;
  outfile << m << " " << m << " " << nb_non_zero << std::endl;

  for (UInt i = 0; i < nb_non_zero; ++i) {
    outfile << irn.values[i] << " " << jcn.values[i] << " 1" << std::endl;
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SparseMatrix::saveMatrix(const std::string & filename) {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.open(filename.c_str());

  outfile << "%%MatrixMarket matrix coordinate real";

  if(sparse_matrix_type == _symmetric) outfile << " symmetric";
  else outfile << " general";
  outfile << std::endl;

  UInt m = size;
  outfile << m << " " << m << " " << nb_non_zero << std::endl;

  for (UInt i = 0; i < nb_non_zero; ++i) {
    outfile << irn.values[i] << " " << jcn.values[i] << " " << a.values[i] << std::endl;
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Vector<Real> & operator*=(Vector<Real> & vect, const SparseMatrix & mat) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT((vect.getSize()*vect.getNbComponent() == mat.getSize()) &&
		      (vect.getNbComponent() == mat.getNbDegreOfFreedom()),
		      "The size of the matrix and the vector do not match");

  UInt nb_non_zero = mat.getNbNonZero();
  Real * tmp = new Real [vect.getNbComponent() * vect.getSize()];
  std::fill_n(tmp, vect.getNbComponent() * vect.getSize(), 0);

  Int * i_val  = mat.getIRN().values;
  Int * j_val  = mat.getJCN().values;
  Real * a_val = mat.getA().values;

  Real * vect_val = vect.values;

  for (UInt k = 0; k < nb_non_zero; ++k) {
    UInt i = *(i_val++);
    UInt j = *(j_val++);
    Real a = *(a_val++);
    tmp[i - 1] += a * vect_val[j - 1];
  }

  memcpy(vect_val, tmp, vect.getNbComponent() * vect.getSize() * sizeof(Real));

  AKANTU_DEBUG_OUT();

  return vect;
}

__END_AKANTU__
