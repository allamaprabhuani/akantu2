/**
 * @file   petsc_matrix.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Jul 21 17:40:41 2014
 *
 * @brief  Implementation of PETSc matrix class
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
#include "petsc_matrix.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
PetscMatrix::PetscMatrix(UInt size,
			   const SparseMatrixType & sparse_matrix_type,
			   UInt nb_degree_of_freedom,
			   const ID & id,
			   const MemoryID & memory_id) :
  PetscMatrix(size, sparse_matrix_type, nb_degree_of_freedom, id, memory_id) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  nb_proc = comm.getNbProc();

  UInt nb_local_nodes = model.getMesh().getNbNodes();

  MatCreate(&comm, &mat);

  MatSetSizes(mat, nb_degree_of_freedom*nb_local_nodes, nb_degree_of_freedom*nb_local_nodes, size, size);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PetscMatrix::PetscMatrix(const PetscMatrix & matrix,
			   const ID & id,
			   const MemoryID & memory_id) :
  PetscMatrix(matrix, id, memory_id) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  nb_proc = comm.getNbProc();

  UInt nb_local_nodes = model.getMesh().getNbNodes();
  MatCreate(&comm, &mat);
  MatSetSizes(mat, nb_degree_of_freedom*nb_local_nodes, nb_degree_of_freedom*nb_local_nodes, size, size);


  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
PetscMatrix::~PetscMatrix() {
  AKANTU_DEBUG_IN();

  MatDestroy(&mat);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PetscMatrix::createAkantuToPetscMap(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int nb_proc = comm.getNbProc();
  Int rank = comm.whoAmI();

  this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

  //initialize vector to store the number of local and master nodes for each processor
  Vector<UInt> master_local_ndofs_per_proc(nb_proc);

  //find the number of dofs corresponding to master or local nodes and store their akantu global dof number
  Real nb_dofs = dof_sychronizer.getNbDOFs();
  Real nb_local_master_dofs = 0;

  Array<Real> local_master_dofs(nb_dofs);
  Array<Real>::scalar_iterator it_dof = local_master_dofs.begin();

  Array<Real> dof_types = dof_sychronizer.getDOFTypes();
  Array<Real>::scalar_iterator it_dof_type = dof_types.begin();

  for (UInt i = 0; i <nb_dofs; ++i; ++it_dof_type) {
    if (*it_dof_type == -2 || *it_dof_type == -1 ) {
      *it_dof = dof_sychronizer.dof_global_ids(i);
      ++it_dof;
      ++nb_local_master_dofs;
    }
  }

  local_master_dofs.resize(nb_local_master_dofs);

  // store the nb of master and local dofs on each processor
  master_local_ndofs_per_proc(rank) = nb_local_master_dofs;
  
  //exchange the information among all processors
  comm.allGather(local_master_ndofs_per_proc.storage(), 1);

  //each processor creates a map for his akantu global dofs to the corresponding petsc global dofs
  UInt petsc_start_index = 0;

  for (UInt i = 0; i < rank; ++i) {
    petsc_start_index +=  local_master_ndofs_per_proc(rank);
  }

  Int * local_master_dofs_ptr = local_master_dofs.storage();
  AOCreateBasic(&comm, nb_local_master_dofs,index_akantu, index_petsc, &ao);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void PetscMatrix::buildProfile(const Mesh & mesh, const DOFSynchronizer & dof_synchronizer) {
  AKANTU_DEBUG_IN();

  // if(irn_jcn_to_k) delete irn_jcn_to_k;
  // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;
  clearProfile();

  this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

  coordinate_list_map::iterator irn_jcn_k_it;

  Int * eq_nb_val = dof_synchronizer.getGlobalDOFEquationNumbers().storage();

  Mesh::type_iterator it  = mesh.firstType(mesh.getSpatialDimension(), _not_ghost, _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType (mesh.getSpatialDimension(), _not_ghost, _ek_not_defined);
  for(; it != end; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

    UInt * conn_val = mesh.getConnectivity(*it, _not_ghost).storage();
    Int * local_eq_nb_val = new Int[nb_degree_of_freedom * nb_nodes_per_element];


    for (UInt e = 0; e < nb_element; ++e) {
      Int * tmp_local_eq_nb_val = local_eq_nb_val;
      for (UInt i = 0; i < nb_nodes_per_element; ++i) {
	UInt n = conn_val[i];
	for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
	  *tmp_local_eq_nb_val++ = eq_nb_val[n * nb_degree_of_freedom + d];
	}
	// memcpy(tmp_local_eq_nb_val, eq_nb_val + n * nb_degree_of_freedom, nb_degree_of_freedom * sizeof(Int));
	// tmp_local_eq_nb_val += nb_degree_of_freedom;
      }

      for (UInt i = 0; i < size_mat; ++i) {
	UInt c_irn = local_eq_nb_val[i];
	if(c_irn < size) {
	  UInt j_start = (sparse_matrix_type == _symmetric) ? i : 0;
	  for (UInt j = j_start; j < size_mat; ++j) {
	    UInt c_jcn = local_eq_nb_val[j];
	    if(c_jcn < size) {
	      MatSetValue(this->mat,1,c_irn, 1, c_jcn, INSERT_VALUES);

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

  /// for pbc @todo correct it for parallel
  if(StaticCommunicator::getStaticCommunicator().getNbProc() == 1) {
    for (UInt i = 0; i < size; ++i) {
      KeyCOO irn_jcn = key(i, i);
      irn_jcn_k_it = irn_jcn_k.find(irn_jcn);
      if(irn_jcn_k_it == irn_jcn_k.end()) {
        irn_jcn_k[irn_jcn] = nb_non_zero;
        irn.push_back(i + 1);
        jcn.push_back(i + 1);
        nb_non_zero++;
      }
    }
  }

  a.resize(nb_non_zero);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PetscMatrix::applyBoundaryNormal(Array<bool> & boundary_normal, Array<Real> & EulerAngles, Array<Real> & rhs, const Array<Real> & matrix, Array<Real> & rhs_rotated) {
    AKANTU_DEBUG_IN();


    const UInt dim = nb_degree_of_freedom;
    const UInt nb_nodes = boundary_normal.getSize();
    Matrix<Real> Ti(dim, dim);
    Matrix<Real> Tj(dim, dim);
    Matrix<Real> small_matrix(dim, dim);
    Matrix<Real> Ti_small_matrix(dim, dim);
    Matrix<Real> Ti_small_matrix_Tj(dim, dim);
    Matrix<Real> small_rhs(dim, 1);
    Matrix<Real> Ti_small_rhs(dim, 1);
    //Transformation matrix based on euler angles X_1Y_2Z_3

    const DOFSynchronizer::GlobalEquationNumberMap & local_eq_num_to_global = dof_synchronizer->getGlobalEquationNumberToLocal();

    Int * irn_val = irn.storage();
    Int * jcn_val = jcn.storage();
    Int * eq_nb_val = dof_synchronizer->getGlobalDOFEquationNumbers().storage();
    Int * eq_nb_rhs_val = dof_synchronizer->getLocalDOFEquationNumbers().storage();
    
    Array<bool> * nodes_rotated = new Array<bool > (nb_nodes, dim, "nodes_rotated");
    nodes_rotated->clear();
    
    Array<Int> * local_eq_nb_val_node_i = new Array<Int > (1, nb_degree_of_freedom, "local_eq_nb_val_node_i");
    local_eq_nb_val_node_i->clear();
    Array<Int> * local_eq_nb_val_node_j = new Array<Int > (1, nb_degree_of_freedom, "local_eq_nb_val_node_j");
    local_eq_nb_val_node_j->clear();

    for (UInt i = 0; i < nb_non_zero; ++i) {
        UInt ni = local_eq_num_to_global.find(*irn_val - 1)->second;
        UInt node_i = (ni - ni % nb_degree_of_freedom) / nb_degree_of_freedom;
        UInt nj = local_eq_num_to_global.find(*jcn_val - 1)->second;
        UInt node_j = (nj - nj % nb_degree_of_freedom) / nb_degree_of_freedom;
        bool constrain_ij = false;
        for (UInt j = 0; j < dim; j++){
            if ( (boundary_normal(node_i, j) || boundary_normal(node_j, j)) ) {
                constrain_ij = true;
                break;
            }
        }
        
        if(constrain_ij){
            if(dim==2){
                Real Theta=EulerAngles(node_i,0);
                Ti(0,0)=cos(Theta);
                Ti(0,1)=-sin(Theta);
                Ti(1,1)=cos(Theta);
                Ti(1,0)=sin(Theta);
                
                Theta=EulerAngles(node_j,0);
                Tj(0,0)=cos(Theta);
                Tj(0,1)=-sin(Theta);
                Tj(1,1)=cos(Theta);
                Tj(1,0)=sin(Theta);
            }
            else if(dim==3){
                Real Theta_x=EulerAngles(node_i,0);
                Real Theta_y=EulerAngles(node_i,1);
                Real Theta_z=EulerAngles(node_i,2);
                
                Ti(0,0)=cos(Theta_y)*cos(Theta_z);
                Ti(0,1)=-cos(Theta_y)*sin(Theta_z);
                Ti(0,2)=sin(Theta_y);
                Ti(1,0)=cos(Theta_x)*sin(Theta_z)+cos(Theta_z)*sin(Theta_x)*sin(Theta_y);
                Ti(1,1)=cos(Theta_x)*cos(Theta_z)-sin(Theta_x)*sin(Theta_y)*sin(Theta_z);
                Ti(1,2)=-cos(Theta_y)*sin(Theta_x);
                Ti(2,0)=sin(Theta_x)*sin(Theta_z)-cos(Theta_x)*cos(Theta_z)*sin(Theta_y);
                Ti(2,1)=cos(Theta_z)*sin(Theta_x)+cos(Theta_x)*sin(Theta_y)*sin(Theta_z);
                Ti(2,2)=cos(Theta_x)*cos(Theta_y);
                
                Theta_x=EulerAngles(node_j,0);
                Theta_y=EulerAngles(node_j,1);
                Theta_z=EulerAngles(node_j,2);
                
                Tj(0,0)=cos(Theta_y)*cos(Theta_z);
                Tj(0,1)=-cos(Theta_y)*sin(Theta_z);
                Tj(0,2)=sin(Theta_y);
                Tj(1,0)=cos(Theta_x)*sin(Theta_z)+cos(Theta_z)*sin(Theta_x)*sin(Theta_y);
                Tj(1,1)=cos(Theta_x)*cos(Theta_z)-sin(Theta_x)*sin(Theta_y)*sin(Theta_z);
                Tj(1,2)=-cos(Theta_y)*sin(Theta_x);
                Tj(2,0)=sin(Theta_x)*sin(Theta_z)-cos(Theta_x)*cos(Theta_z)*sin(Theta_y);
                Tj(2,1)=cos(Theta_z)*sin(Theta_x)+cos(Theta_x)*sin(Theta_y)*sin(Theta_z);
                Tj(2,2)=cos(Theta_x)*cos(Theta_y);
            }
            for (UInt j = 0; j < nb_degree_of_freedom; ++j){
                local_eq_nb_val_node_i->storage()[j] = eq_nb_val[node_i * nb_degree_of_freedom + j];
                local_eq_nb_val_node_j->storage()[j] = eq_nb_val[node_j * nb_degree_of_freedom + j];
            }
            for (UInt j = 0; j < nb_degree_of_freedom; ++j) {
                for (UInt k = 0; k < nb_degree_of_freedom; ++k) {
                    KeyCOO jcn_irn = key(local_eq_nb_val_node_i->storage()[j], local_eq_nb_val_node_j->storage()[k]);
                    coordinate_list_map::iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);
                    small_matrix(j, k) = matrix.storage()[irn_jcn_k_it->second];
                }
                small_rhs(j,0)=rhs.storage()[eq_nb_rhs_val[node_i*dim+j]];
            }
            
            Ti_small_matrix.mul < false, false > (Ti, small_matrix);
            Ti_small_matrix_Tj.mul < false, true> (Ti_small_matrix, Tj);
            Ti_small_rhs.mul<false, false>(Ti, small_rhs);
            
            /*for (UInt j = 0; j < nb_degree_of_freedom; ++j) {
                for (UInt k = 0; k < nb_degree_of_freedom; ++k) {
                    KeyCOO jcn_irn = key(local_eq_nb_val_node_i->storage()[j], local_eq_nb_val_node_j->storage()[k]);
                    coordinate_list_map::iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);
                    a.storage()[irn_jcn_k_it->second] = T_small_matrix_T(j,k);
                }
                if(node_constrain==node_i && !( nodes_rotated->storage()[node_i]))
                    rhs.storage()[eq_nb_rhs_val[node_i*dim+j]]=T_small_rhs(j,0);
            }*/
            KeyCOO jcn_irn = key(ni, nj);
            coordinate_list_map::iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);
            a.storage()[irn_jcn_k_it->second] = Ti_small_matrix_Tj(ni % nb_degree_of_freedom, nj % nb_degree_of_freedom);
            //this->saveMatrix("stiffness_partial.out");
            if (!(nodes_rotated->storage()[eq_nb_rhs_val[node_i * dim + ni % nb_degree_of_freedom]])) {
                rhs_rotated.storage()[eq_nb_rhs_val[node_i * dim + ni % nb_degree_of_freedom]] = Ti_small_rhs(ni % nb_degree_of_freedom, 0);
                nodes_rotated->storage()[eq_nb_rhs_val[node_i * dim + ni % nb_degree_of_freedom]] = true;
            }
        }
        else{
            if (!(nodes_rotated->storage()[eq_nb_rhs_val[node_i * dim + ni % nb_degree_of_freedom]])) {
                rhs_rotated.storage()[eq_nb_rhs_val[node_i * dim + ni % nb_degree_of_freedom]] = rhs.storage()[eq_nb_rhs_val[node_i * dim + ni % nb_degree_of_freedom]];
                nodes_rotated->storage()[eq_nb_rhs_val[node_i * dim + ni % nb_degree_of_freedom]] = true;
            }
        }
        irn_val++;
        jcn_val++;
    }

    delete local_eq_nb_val_node_i;
    delete local_eq_nb_val_node_j;
    delete nodes_rotated;
    this->saveMatrix("stiffness_rotated.out");
    applyBoundary(boundary_normal);

    /*Real norm = 0;
    UInt n_zeros = 0;
    for (UInt j = 0; j < nb_degree_of_freedom; ++j) {
        norm += Normal.storage()[j] * Normal.storage()[j];
        if (std::abs(Normal.storage()[j]) < Math::getTolerance())
            n_zeros++;
    }
    norm = sqrt(norm);
    AKANTU_DEBUG_ASSERT(norm > Math::getTolerance(), "The normal is not right");
    for (UInt j = 0; j < nb_degree_of_freedom; ++j)
        Normal.storage()[j] /= norm;
    
    if (n_zeros == nb_degree_of_freedom - 1) {
        UInt nb_nodes = boundary_normal.getSize();
        for (UInt i = 0; i < nb_nodes; ++i) {
            if (boundary_normal(i, 0))
                for (UInt j = 0; j < nb_degree_of_freedom; ++j)
                    boundary(i, j) += Normal(0, j);
        }
    } else {

        const DOFSynchronizer::GlobalEquationNumberMap & local_eq_num_to_global = dof_synchronizer->getGlobalEquationNumberToLocal();

        Int * irn_val = irn.storage();
        Int * eq_nb_val = dof_synchronizer->getGlobalDOFEquationNumbers().storage();
        Array<Int> * local_eq_nb_val = new Array<Int > (1, nb_degree_of_freedom, "local_eq_nb_val");
        local_eq_nb_val->clear();
        UInt nb_nodes = boundary_normal.getSize();
        Array<bool> * node_set = new Array<bool > (nb_nodes, 1, "node_set");
        node_set->clear();

        for (UInt i = 0; i < nb_non_zero; ++i) {
            UInt ni = local_eq_num_to_global.find(*irn_val - 1)->second;
            UInt node_i = (ni - ni % nb_degree_of_freedom) / nb_degree_of_freedom;
            if (boundary_normal.storage()[ni]) {
                if ((!number.storage()[ni]) && (!node_set->storage()[node_i])) {
                    UInt normal_component_i = ni % nb_degree_of_freedom; //DOF of node node_i to be removed
                    if (std::abs(Normal.storage()[normal_component_i]) > Math::getTolerance()) {
                        for (UInt j = 0; j < nb_degree_of_freedom; ++j)
                            local_eq_nb_val->storage()[j] = eq_nb_val[node_i * nb_degree_of_freedom + j];

                        UInt DOF_remove = local_eq_nb_val->storage()[normal_component_i];
                        KeyCOO jcn_irn = key(DOF_remove, DOF_remove);
                        coordinate_list_map::iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);

                        Real aii = a.storage()[irn_jcn_k_it->second];
                        Real normal_constant = (1 - aii) / (Normal.storage()[normal_component_i] * Normal.storage()[normal_component_i]);

                        for (UInt j = 0; j < nb_degree_of_freedom; ++j) {
                            UInt line_j = local_eq_nb_val->storage()[j];
                            for (UInt k = 0; k < nb_degree_of_freedom; ++k) {
                                UInt column_k = local_eq_nb_val->storage()[k];
                                jcn_irn = key(line_j, column_k);
                                if ((line_j == column_k) && (line_j == DOF_remove))
                                    a.storage()[irn_jcn_k_it->second] = 1;
                                else
                                    a.storage()[irn_jcn_k_it->second] += Normal.storage()[j] * Normal.storage()[k] * normal_constant;
                            }
                        }

                        number.storage()[ni]++;
                        node_set->storage()[node_i] = false;
                    }
                }
            }
            irn_val++;
        }

        delete local_eq_nb_val;
        delete node_set;
    }*/
    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PetscMatrix::applyBoundary(const Array<bool> & boundary, Real block_val) {
  AKANTU_DEBUG_IN();

  const DOFSynchronizer::GlobalEquationNumberMap & local_eq_num_to_global = dof_synchronizer->getGlobalEquationNumberToLocal();
  Int * irn_val = irn.storage();
  Int * jcn_val = jcn.storage();
  Real * a_val   = a.storage();

  for (UInt i = 0; i < nb_non_zero; ++i) {
    UInt ni = local_eq_num_to_global.find(*irn_val - 1)->second;
    UInt nj = local_eq_num_to_global.find(*jcn_val - 1)->second;
    if(boundary.storage()[ni]  || boundary.storage()[nj]) {
      if (*irn_val != *jcn_val) {
        *a_val = 0;
      } else {
        if(dof_synchronizer->getDOFTypes()(ni) >= 0) {
          *a_val = 0;
        } else {
          *a_val = block_val;
        }
      }
    }
    irn_val++; jcn_val++; a_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PetscMatrix::removeBoundary(const Array<bool> & boundary) {
  AKANTU_DEBUG_IN();

  if(irn_save) delete irn_save;
  if(jcn_save) delete jcn_save;

  irn_save = new Array<Int>(irn, true);
  jcn_save = new Array<Int>(jcn, true);

  UInt n = boundary.getSize()*boundary.getNbComponent();

  UInt * perm = new UInt[n];

  size_save = size;
  size = 0;
  for (UInt i = 0; i < n; ++i) {
    if(!boundary.storage()[i]) {
      perm[i] = size;
      //      std::cout <<  "perm["<< i <<"] = " << size << std::endl;
      size++;
    }
  }

  for (UInt i = 0; i < nb_non_zero;) {
    if(boundary.storage()[irn(i) - 1] || boundary.storage()[jcn(i) - 1]) {
      irn.erase(i);
      jcn.erase(i);
      a.erase(i);
      nb_non_zero--;
    } else {
      irn(i) = perm[irn(i) - 1] + 1;
      jcn(i) = perm[jcn(i) - 1] + 1;
      i++;
    }
  }

  delete [] perm;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PetscMatrix::restoreProfile() {
  AKANTU_DEBUG_IN();

  irn.resize(irn_save->getSize());
  jcn.resize(jcn_save->getSize());

  nb_non_zero = irn.getSize();
  a.resize(nb_non_zero);
  size = size_save;

  memcpy(irn.storage(), irn_save->storage(), irn.getSize()*sizeof(Int));
  memcpy(jcn.storage(), jcn_save->storage(), jcn.getSize()*sizeof(Int));

  delete irn_save; irn_save = NULL;
  delete jcn_save; jcn_save = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PetscMatrix::saveProfile(const std::string & filename) const {
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
    outfile << irn.storage()[i] << " " << jcn.storage()[i] << " 1" << std::endl;
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PetscMatrix::saveMatrix(const std::string & filename) const {
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
Array<Real> & operator*=(Array<Real> & vect, const PetscMatrix & mat) {
  AKANTU_DEBUG_IN();

  Array<Real> result 
  MatMult(this->mat, vect, )

  // AKANTU_DEBUG_ASSERT((vect.getSize()*vect.getNbComponent() == mat.getSize()) &&
  // 		      (vect.getNbComponent() == mat.getNbDegreOfFreedom()),
  // 		      "The size of the matrix and the vector do not match");

  const PetscMatrixType & sparse_matrix_type = mat.getPetscMatrixType();
  DOFSynchronizer * dof_synchronizer = mat.getDOFSynchronizerPointer();

  UInt nb_non_zero = mat.getNbNonZero();
  Real * tmp = new Real [vect.getNbComponent() * vect.getSize()];
  std::fill_n(tmp, vect.getNbComponent() * vect.getSize(), 0);

  Int * i_val  = mat.getIRN().storage();
  Int * j_val  = mat.getJCN().storage();
  Real * a_val = mat.getA().storage();

  Real * vect_val = vect.storage();

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
    if((sparse_matrix_type == _symmetric) && (local_i != local_j))
      tmp[local_j] += a * vect_val[local_i];
  }

  memcpy(vect_val, tmp, vect.getNbComponent() * vect.getSize() * sizeof(Real));
  delete [] tmp;

  if(dof_synchronizer)
    dof_synchronizer->reduceSynchronize<AddOperation>(vect);

  AKANTU_DEBUG_OUT();

  return vect;
}

/* -------------------------------------------------------------------------- */
void PetscMatrix::copyContent(const PetscMatrix & matrix) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
		      "The two matrices don't have the same profiles");

  MatCopy(this->mat, matrix->mat, SAME_NONZERO_PATTERN);

  AKANTU_DEBUG_OUT();
}

///* -------------------------------------------------------------------------- */
//void PetscMatrix::copyProfile(const PetscMatrix & matrix) {
//  AKANTU_DEBUG_IN();
//  irn = matrix.irn;
//  jcn = matrix.jcn;
//  nb_non_zero = matrix.nb_non_zero;
//  irn_jcn_k = matrix.irn_jcn_k;
//  a.resize(nb_non_zero);
//  AKANTU_DEBUG_OUT();
//}

/* -------------------------------------------------------------------------- */
void PetscMatrix::add(const SparseMatrix & matrix, Real alpha) {
  AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
		      "The two matrices don't have the same profile");

  // get values and local indices of matrix that has to be added
  UInt i_global;
  UInt j_global;

  Array<UInt>::scalar_iterator i_global_it = i_global.begin();
  Array<UInt>::scalar_iterator i_global_end = i_global.end();
  Array<UInt>::scalar_iterator j_global_it = j_global.begin();

  Int * i = matrix.irn.storage();
  Int *j = matrix.jcn.storage();
  Real * a_val = matrix.a.storage();

  for (; i != i_global_end; ++i_global_it; ++j_global_it; ++i; ++j; ++a_val) {
    // get the akantu global dof index
    i_global = dof_synchronizer->getDOFGlobalID(*i);
    j_global = dof_synchronizer->getDOFGlobalID(*j);

    // get the petsc gloobal dof index
    AOApplicationToPetsc(&ao,2,{i_global,j_global})
    MatSetValue(*mat,i,j,*a_val,INSERT_VALUES); 
  }

 
}


/* -------------------------------------------------------------------------- */
void PetscMatrix::lump(Array<Real> & lumped) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT((lumped.getNbComponent() == nb_degree_of_freedom),
		      "The size of the matrix and the vector do not match");

  UInt vect_size = size / nb_degree_of_freedom;
  if(dof_synchronizer) vect_size = dof_synchronizer->getNbDOFs() / nb_degree_of_freedom;

  lumped.resize(vect_size);
  lumped.clear();

  Int * i_val  = irn.storage();
  Int * j_val  = jcn.storage();
  Real * a_val = a.storage();

  Real * vect_val = lumped.storage();

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
    dof_synchronizer->reduceSynchronize<AddOperation>(lumped);

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
