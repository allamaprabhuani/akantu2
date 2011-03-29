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

// /* -------------------------------------------------------------------------- */
// SparseMatrix::SparseMatrix(const Mesh & mesh,
// 			   const SparseMatrixType & sparse_matrix_type,
// 			   UInt nb_degre_of_freedom,
// 			   const SparseMatrixID & id,
// 			   const MemoryID & memory_id) :
//   Memory(memory_id), id(id),
//   sparse_matrix_type(sparse_matrix_type),
//   nb_degre_of_freedom(nb_degre_of_freedom),
//   mesh(&mesh),
//   nb_non_zero(0),
//   irn(0,1,"irn"), jcn(0,1,"jcn"), a(0,1,"A"),
//   irn_jcn_to_k(NULL) {
//   AKANTU_DEBUG_IN();

//   size = mesh.getNbGlobalNodes()*nb_degre_of_freedom;

//   for(UInt t = _not_defined; t < _max_element_type; ++t) {
//     this->element_to_sparse_profile[t] = NULL;
//   }

//   StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
//   nb_proc = comm->getNbProc();

//   irn_save = NULL;
//   jcn_save = NULL;

//   irn_jcn_k = new coordinate_list_map;

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(UInt size,
			   const SparseMatrixType & sparse_matrix_type,
			   UInt nb_degre_of_freedom,
			   const SparseMatrixID & id,
			   const MemoryID & memory_id) :
  Memory(memory_id), id(id),
  sparse_matrix_type(sparse_matrix_type),
  nb_degre_of_freedom(nb_degre_of_freedom),
  size(size),
  nb_non_zero(0),
  irn(0,1,"irn"), jcn(0,1,"jcn"), a(0,1,"A") {
  AKANTU_DEBUG_IN();

  // for(UInt t = _not_defined; t < _max_element_type; ++t) {
  //   this->element_to_sparse_profile[t] = NULL;
  // }

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
  size(matrix.getSize()), nb_non_zero(matrix.getNbNonZero()),
  irn(matrix.getIRN(), true), jcn(matrix.getJCN(), true), a(matrix.getA(), true) {
  AKANTU_DEBUG_IN();

  irn_save = NULL;
  jcn_save = NULL;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
SparseMatrix::~SparseMatrix() {
  AKANTU_DEBUG_IN();

  //  if (irn_jcn_to_k) delete irn_jcn_to_k;

  if(irn_save) delete irn_save;
  if(jcn_save) delete jcn_save;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::buildProfile(const Mesh & mesh, const Vector<Int> & equation_number) {
  AKANTU_DEBUG_IN();

  // if(irn_jcn_to_k) delete irn_jcn_to_k;
  // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;
  clearProfile();

  coordinate_list_map::iterator irn_jcn_k_it;

  // UInt * global_nodes_ids_val = NULL;
  // if(nb_proc > 1) global_nodes_ids_val = mesh.getGlobalNodesIds().values;

  Int * eq_nb_val = equation_number.values;

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if (Mesh::getSpatialDimension(*it) != mesh.getSpatialDimension()) continue;

    UInt nb_element = mesh.getNbElement(*it);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt size_mat = nb_nodes_per_element * nb_degre_of_freedom;

    UInt * conn_val = mesh.getConnectivity(*it).values;
    Int * local_eq_nb_val = new Int[nb_degre_of_freedom * nb_nodes_per_element];


    for (UInt e = 0; e < nb_element; ++e) {
      Int * tmp_local_eq_nb_val = local_eq_nb_val;
      for (UInt i = 0; i < nb_nodes_per_element; ++i) {
	UInt n = conn_val[i];
	memcpy(tmp_local_eq_nb_val, eq_nb_val + n * nb_degre_of_freedom, nb_degre_of_freedom * sizeof(Int));
	tmp_local_eq_nb_val += nb_degre_of_freedom;
      }

      for (UInt i = 0; i < size_mat; ++i) {
	UInt c_irn = local_eq_nb_val[i];
	if(c_irn < size) {
	  UInt j_start = (sparse_matrix_type == _symmetric) ? i : 0;
	  for (UInt j = j_start; j < size_mat; ++j) {
	    UInt c_jcn = local_eq_nb_val[j];
	    if(c_jcn < size) {
	      KeyCOO irn_jcn = key(c_irn, c_jcn);
	      irn_jcn_k_it = irn_jcn_k.find(irn_jcn);

	      if (irn_jcn_k_it == irn_jcn_k.end()) {
		irn_jcn_k[irn_jcn] = nb_non_zero;
		irn.push_back(c_irn + 1);
		jcn.push_back(c_jcn + 1);
		nb_non_zero++;
	      }
	    }
	  }
	}
      }
      conn_val += nb_nodes_per_element;
    }

    delete [] local_eq_nb_val;
    // for (UInt e = 0; e < nb_element; ++e) {                                     // loop on element
    //   for (UInt j = 0; j < nb_nodes_per_element; ++j) {                         // loop on local column
    // 	UInt n_j = (nb_proc == 1) ? conn_val[j] : global_nodes_ids_val[conn_val[j]];
    // 	UInt c_jcn = n_j * nb_degre_of_freedom;

    // 	for (UInt d_j = 0; d_j < nb_degre_of_freedom; ++d_j, ++c_jcn) {         // loop on degre of freedom
    // 	  for (UInt i = 0; i < nb_nodes_per_element; ++i) {                                    // loop on rows
    // 	    UInt n_i = (nb_proc == 1) ? conn_val[i] : global_nodes_ids_val[conn_val[i]];
    // 	    UInt c_irn = n_i * nb_degre_of_freedom;

    // 	    for (UInt d_i = 0; d_i < nb_degre_of_freedom; ++d_i, ++c_irn) {     // loop on degre of freedom
    // 	      UInt irn_jcn;
    // 	      irn_jcn = key(c_irn, c_jcn);

    // 	      irn_jcn_k_it = irn_jcn_k.find(irn_jcn);

    // 	      if (irn_jcn_k_it == irn_jcn_k.end()) {
    // 		irn_jcn_k[irn_jcn] = nb_non_zero;
    // 		irn.push_back(c_irn + 1);
    // 		jcn.push_back(c_jcn + 1);
    // 		nb_non_zero++;
    // 	      }
    // 	    }
    // 	  }
    // 	}
    //   }
    //   conn_val += nb_nodes_per_element;
    // }
  }

  a.resize(nb_non_zero);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::applyBoundary(const Vector<bool> & boundary,
				 const unordered_map<UInt, UInt>::type & local_eq_num_to_global) {
  AKANTU_DEBUG_IN();

  Int * irn_val = irn.values;
  Int * jcn_val = jcn.values;
  Real * a_val   = a.values;

  for (UInt i = 0; i < nb_non_zero; ++i) {
    UInt ni = local_eq_num_to_global.find(*irn_val - 1)->second;
    UInt nj = local_eq_num_to_global.find(*jcn_val - 1)->second;
    if(boundary.values[ni]  || boundary.values[nj]) {
     if (*irn_val != *jcn_val) *a_val = 0;
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




  outfile.open("toto.mtx");

  outfile << "%%MatrixMarket matrix coordinate pattern";

  if(sparse_matrix_type == _symmetric) outfile << " symmetric";
  else outfile << " general";
  outfile << std::endl;

  outfile << m << " " << m << " " << nb_non_zero << std::endl;

  coordinate_list_map::const_iterator it;
  for (it = irn_jcn_k.begin(); it != irn_jcn_k.end(); ++it) {
    //    outfile << it->first / size + 1 << " " << it->first % size + 1 << " 1" << std::endl;
    outfile << it->first.first + 1 << " " << it->first.second + 1<< " 1" << std::endl;
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SparseMatrix::saveMatrix(const std::string & filename) {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.precision(std::numeric_limits<Real>::digits10);

  outfile.open(filename.c_str());

  outfile << "%%MatrixMarket matrix coordinate real";

  if(sparse_matrix_type == _symmetric) outfile << " symmetric";
  else outfile << " general";
  outfile << std::endl;

  outfile << size << " " << size << " " << nb_non_zero << std::endl;

  coordinate_list_map::const_iterator it;
  for (it = irn_jcn_k.begin(); it != irn_jcn_k.end(); ++it) {
    outfile << it->first.first + 1 << " " << it->first.second + 1 << " " << a.values[it->second] << std::endl;
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
