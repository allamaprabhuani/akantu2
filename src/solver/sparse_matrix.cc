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
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(UInt size,
			   const SparseMatrixType & sparse_matrix_type,
			   UInt nb_degre_of_freedom,
			   const ID & id,
			   const MemoryID & memory_id) :
  Memory(memory_id), id(id),
  sparse_matrix_type(sparse_matrix_type),
  nb_degre_of_freedom(nb_degre_of_freedom),
  size(size),
  nb_non_zero(0),
  irn(0,1,"irn"), jcn(0,1,"jcn"), a(0,1,"A") {
  AKANTU_DEBUG_IN();

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  nb_proc = comm->getNbProc();
  dof_synchronizer = NULL;

  irn_save = NULL;
  jcn_save = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(const SparseMatrix & matrix,
			   const ID & id,
			   const MemoryID & memory_id) :
  Memory(memory_id), id(id),
  sparse_matrix_type(matrix.sparse_matrix_type),
  nb_degre_of_freedom(matrix.nb_degre_of_freedom),
  size(matrix.size),
  nb_proc(matrix.nb_proc),
  nb_non_zero(matrix.nb_non_zero),
  irn(matrix.irn, true), jcn(matrix.jcn, true), a(matrix.getA(), true),
  irn_jcn_k(matrix.irn_jcn_k) {
  AKANTU_DEBUG_IN();

  size_save = 0;
  irn_save = NULL;
  jcn_save = NULL;
  dof_synchronizer = matrix.dof_synchronizer;

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
void SparseMatrix::buildProfile(const Mesh & mesh, const DOFSynchronizer & dof_synchronizer) {
  AKANTU_DEBUG_IN();

  // if(irn_jcn_to_k) delete irn_jcn_to_k;
  // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;
  clearProfile();

  this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

  coordinate_list_map::iterator irn_jcn_k_it;

  Int * eq_nb_val = dof_synchronizer.getGlobalDOFEquationNumbers().values;

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if (Mesh::getSpatialDimension(*it) != mesh.getSpatialDimension()) continue;

    UInt nb_element = mesh.getNbElement(*it);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt size_mat = nb_nodes_per_element * nb_degre_of_freedom;

    UInt * conn_val = mesh.getConnectivity(*it, _not_ghost).values;
    Int * local_eq_nb_val = new Int[nb_degre_of_freedom * nb_nodes_per_element];


    for (UInt e = 0; e < nb_element; ++e) {
      Int * tmp_local_eq_nb_val = local_eq_nb_val;
      for (UInt i = 0; i < nb_nodes_per_element; ++i) {
	UInt n = conn_val[i];
	for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	  *tmp_local_eq_nb_val++ = eq_nb_val[n * nb_degre_of_freedom + d];
	}
	// memcpy(tmp_local_eq_nb_val, eq_nb_val + n * nb_degre_of_freedom, nb_degre_of_freedom * sizeof(Int));
	// tmp_local_eq_nb_val += nb_degre_of_freedom;
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
  }

  for (UInt i = 0; i < size; ++i) {
    KeyCOO irn_jcn = key(i, i);
    coordinate_list_map::const_iterator irn_jcn_k_it = irn_jcn_k.find(irn_jcn);
    if(irn_jcn_k_it == irn_jcn_k.end()) {
      irn_jcn_k[irn_jcn] = nb_non_zero;
      irn.push_back(i + 1);
      jcn.push_back(i + 1);
      nb_non_zero++;
    }
  }

  a.resize(nb_non_zero);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::applyBoundary(const Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  const DOFSynchronizer::GlobalEquationNumberMap & local_eq_num_to_global = dof_synchronizer->getGlobalEquationNumberToLocal();
  Int * irn_val = irn.values;
  Int * jcn_val = jcn.values;
  Real * a_val   = a.values;

  for (UInt i = 0; i < nb_non_zero; ++i) {
    UInt ni = local_eq_num_to_global.find(*irn_val - 1)->second;
    UInt nj = local_eq_num_to_global.find(*jcn_val - 1)->second;
    if(boundary.values[ni]  || boundary.values[nj]) {
     if (*irn_val != *jcn_val) *a_val = 0;
     else {
       if(dof_synchronizer->getDOFTypes()(ni) >= 0) *a_val = 0;
       else *a_val = 1;
     }
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
      //      std::cout <<  "perm["<< i <<"] = " << size << std::endl;
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
void SparseMatrix::saveProfile(const std::string & filename) const {
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
void SparseMatrix::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.precision(std::numeric_limits<Real>::digits10);

  outfile.open(filename.c_str());

  outfile << "%%MatrixMarket matrix coordinate real";

  if(sparse_matrix_type == _symmetric) outfile << " symmetric";
  else outfile << " general";
  outfile << std::endl;

  outfile << size << " " << size << " " << nb_non_zero << std::endl;

  for (UInt i = 0; i < nb_non_zero; ++i) {
    outfile << irn(i) << " " << jcn(i) << " " << a(i) << std::endl;
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Vector<Real> & operator*=(Vector<Real> & vect, const SparseMatrix & mat) {
  AKANTU_DEBUG_IN();

  // AKANTU_DEBUG_ASSERT((vect.getSize()*vect.getNbComponent() == mat.getSize()) &&
  // 		      (vect.getNbComponent() == mat.getNbDegreOfFreedom()),
  // 		      "The size of the matrix and the vector do not match");

  const SparseMatrixType & sparse_matrix_type = mat.getSparseMatrixType();
  DOFSynchronizer * dof_synchronizer = mat.getDOFSynchronizerPointer();

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

    UInt local_i = i - 1;
    UInt local_j = j - 1;
    if(dof_synchronizer) {
      local_i = dof_synchronizer->getDOFLocalID(local_i);
      local_j = dof_synchronizer->getDOFLocalID(local_j);
    }

    tmp[local_i] += a * vect_val[local_j];
    if(sparse_matrix_type == _symmetric && (local_i != local_j))
      tmp[local_j] += a * vect_val[local_i];
  }

  memcpy(vect_val, tmp, vect.getNbComponent() * vect.getSize() * sizeof(Real));
  delete [] tmp;

  if(dof_synchronizer)
    dof_synchronizer->reduceSynchronize<AddOperation<Real> >(vect);

  AKANTU_DEBUG_OUT();

  return vect;
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::copyContent(const SparseMatrix & matrix) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
		      "The to matrix don't have the same profiles");
  memcpy(a.values, matrix.getA().values, nb_non_zero * sizeof(Real));
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrix::add(const SparseMatrix & matrix, Real alpha) {
  AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
		      "The to matrix don't have the same profiles");

  Real * a_val = a.values;
  Real * b_val = matrix.a.values;

  for (UInt n = 0; n < nb_non_zero; ++n) {
    *a_val++ += alpha * *b_val++;
  }
}


/* -------------------------------------------------------------------------- */
void SparseMatrix::lump(Vector<Real> & lumped) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT((lumped.getNbComponent() == nb_degre_of_freedom),
		      "The size of the matrix and the vector do not match");

  UInt vect_size = size / nb_degre_of_freedom;
  if(dof_synchronizer) vect_size = dof_synchronizer->getNbDOFs() / nb_degre_of_freedom;

  lumped.resize(vect_size);
  lumped.clear();

  Int * i_val  = irn.values;
  Int * j_val  = jcn.values;
  Real * a_val = a.values;

  Real * vect_val = lumped.values;

  for (UInt k = 0; k < nb_non_zero; ++k) {
    UInt i = *(i_val++);
    UInt j = *(j_val++);
    Real a = *(a_val++);

    UInt local_i = i - 1;
    UInt local_j = j - 1;
    if(dof_synchronizer) {
      local_i = dof_synchronizer->getDOFLocalID(local_i);
      local_j = dof_synchronizer->getDOFLocalID(local_j);
    }

    vect_val[local_i] += a;
    if(sparse_matrix_type == _symmetric && (i != j))
      vect_val[local_j] += a;
  }

  if(dof_synchronizer)
    dof_synchronizer->reduceSynchronize<AddOperation<Real> >(lumped);

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
