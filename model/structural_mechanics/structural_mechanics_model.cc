/**
 * @file   structural_mechanics_model.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu May  5 15:52:38 2011
 *
 * @brief
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
#include "structural_mechanics_model.hh"
#include "aka_math.hh"
#include "integration_scheme_2nd_order.hh"

#include "static_communicator.hh"
#include "sparse_matrix.hh"
#include "solver.hh"

#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::StructuralMechanicsModel(Mesh & mesh,
						   UInt dim,
						   const ModelID & id,
						   const MemoryID & memory_id) :
  Model(id, memory_id),
  stiffness_matrix(NULL),
  solver(NULL),
  spatial_dimension(dim) {
  AKANTU_DEBUG_IN();

  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();
  registerFEMObject<MyFEMType>("SolidMechanicsFEM", mesh, spatial_dimension);

  this->displacement_rotation = NULL;
  this->force_momentum        = NULL;
  this->residual     = NULL;
  this->boundary     = NULL;
  this->increment    = NULL;

  if(spatial_dimension == 2)
    nb_degre_of_freedom = 3;
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::~StructuralMechanicsModel() {
  AKANTU_DEBUG_IN();

  if(solver) delete solver;
  if(stiffness_matrix) delete stiffness_matrix;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initVectors() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEM().getMesh().getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":boundary";
  std::stringstream sstr_incr; sstr_incr << id << ":increment";

  displacement_rotation = &(alloc<Real>(sstr_disp.str(), nb_nodes, nb_degre_of_freedom, REAL_INIT_VALUE));
  force_momentum        = &(alloc<Real>(sstr_forc.str(), nb_nodes, nb_degre_of_freedom, REAL_INIT_VALUE));
  residual     = &(alloc<Real>(sstr_resi.str(), nb_nodes, nb_degre_of_freedom, REAL_INIT_VALUE));
  boundary     = &(alloc<bool>(sstr_boun.str(), nb_nodes, nb_degre_of_freedom, false));
  increment     = &(alloc<Real>(sstr_incr.str(), nb_nodes, nb_degre_of_freedom, REAL_INIT_VALUE));

  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    UInt nb_element = getFEM().getMesh().getNbElement(*it);
    UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(*it);

    if(!element_material.exists(*it, _not_ghost)) {
      std::stringstream sstr_elma; sstr_elma << id << ":element_material:" << *it;
      element_material(*it, _not_ghost) = &(alloc<UInt>(sstr_elma.str(), nb_element, 1, 0));
    }
    if(!stress.exists(*it, _not_ghost)) {
      UInt size = getTangentStiffnessVoigtSize(*it);
      std::stringstream sstr_stress; sstr_stress << id << ":stress:" << *it;
      stress(*it, _not_ghost) = &(alloc<Real>(sstr_stress.str(), nb_element * nb_quadrature_points, size , 0));
    }
  }

  dof_synchronizer = new DOFSynchronizer(getFEM().getMesh(), nb_degre_of_freedom);
  dof_synchronizer->initLocalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initModel() {
  getFEM().initShapeFunctions(_not_ghost);
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initImplicitSolver() {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = getFEM().getMesh();
  UInt nb_nodes = mesh.getNbNodes();

  std::stringstream sstr; sstr << id << ":stiffness_matrix";
  stiffness_matrix = new SparseMatrix(mesh.getNbGlobalNodes() * nb_degre_of_freedom, _symmetric,
				      nb_degre_of_freedom, sstr.str(), memory_id);

  dof_synchronizer->initGlobalDOFEquationNumbers();

  stiffness_matrix->buildProfile(mesh, *dof_synchronizer);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << id << ":solver_stiffness_matrix";
  solver = new SolverMumps(*stiffness_matrix, sstr_solv.str());
  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS
  solver->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt StructuralMechanicsModel::getTangentStiffnessVoigtSize(const ElementType & type) {
  UInt size;
#define GET_TANGENT_STIFFNESS_VOIGT_SIZE(type)	\
  size = getTangentStiffnessVoigtSize<type>();

  AKANTU_BOOST_ELEMENT_SWITCH(GET_TANGENT_STIFFNESS_VOIGT_SIZE);
#undef GET_TANGENT_STIFFNESS_VOIGT_SIZE

  return size;
}
/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  stiffness_matrix->clear();

  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;

#define ASSEMBLE_STIFFNESS_MATRIX(type)		\
    assembleStiffnessMatrix<type>();

    AKANTU_BOOST_ELEMENT_SWITCH(ASSEMBLE_STIFFNESS_MATRIX);
#undef ASSEMBLE_STIFFNESS_MATRIX
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeStressOnQuad() {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;

#define COMPUTE_STRESS_ON_QUAD(type)		\
    computeStressOnQuad<type>();

    AKANTU_BOOST_ELEMENT_SWITCH(COMPUTE_STRESS_ON_QUAD);
#undef COMPUTE_STRESS_ON_QUAD
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::updateResidual() {
 AKANTU_DEBUG_IN();
 residual->copy(*force_momentum);

 Vector<Real> ku(*displacement_rotation, true);
 ku *= *stiffness_matrix;
 *residual -= ku;

 AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::solve() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving an implicit step.");

  UInt nb_nodes = displacement_rotation->getSize();

  /// todo residual = force - Kxr * d_bloq
  stiffness_matrix->applyBoundary(*boundary);

  solver->setRHS(*residual);

  solver->solve(*increment);

  Real * increment_val     = increment->values;
  Real * displacement_val  = displacement_rotation->values;
  bool * boundary_val      = boundary->values;

  for (UInt n = 0; n < nb_nodes * nb_degre_of_freedom; ++n) {
    if(!(*boundary_val)) {
      *displacement_val += *increment_val;
    }

    displacement_val++;
    boundary_val++;
    increment_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool StructuralMechanicsModel::testConvergenceIncrement(Real tolerance) {
  Real error;
  bool tmp = testConvergenceIncrement(tolerance, error);

  AKANTU_DEBUG_INFO("Norm of increment : " << error);

  return tmp;
}

/* -------------------------------------------------------------------------- */
bool StructuralMechanicsModel::testConvergenceIncrement(Real tolerance, Real & error) {
  AKANTU_DEBUG_IN();
  Mesh & mesh= getFEM().getMesh();
  UInt nb_nodes = displacement_rotation->getSize();
  UInt nb_degre_of_freedom = displacement_rotation->getNbComponent();

  Real norm = 0;
  Real * increment_val     = increment->values;
  bool * boundary_val      = boundary->values;

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
      if(!(*boundary_val) && is_local_node) {
	norm += *increment_val * *increment_val;
      }
      boundary_val++;
      increment_val++;
    }
  }

  StaticCommunicator::getStaticCommunicator()->allReduce(&norm, 1, _so_sum);

  error = sqrt(norm);
  AKANTU_DEBUG_ASSERT(!isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (error < tolerance);
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeTangentStiffness<_bernoulli_beam_2>(Vector<Real> & tangent_stiffness_matrix) {
  UInt nb_element                 = getFEM().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(_bernoulli_beam_2);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(_bernoulli_beam_2);


  UInt tangent_size = 2;

  tangent_stiffness_matrix.clear();
  Vector<Real>::iterator<types::Matrix> D = tangent_stiffness_matrix.begin(tangent_size, tangent_size);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = (*element_material(_bernoulli_beam_2, _not_ghost))(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real I = materials[mat].I;
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      (*D)(0,0) = E * A;
      (*D)(1,1) = E * I;
      ++D;
    }
  }
}
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<_bernoulli_beam_2>(Vector<Real> & b) {
  MyFEMType & fem = getFEMClass<MyFEMType>();
  const Vector<Real> & Np  = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 0);
  const Vector<Real> & Mpp = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 1);
  const Vector<Real> & Lpp = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 2);

  UInt nb_element                 = getFEM().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(_bernoulli_beam_2);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(_bernoulli_beam_2);

  Real * Np_val  = Np.values;
  Real * Mpp_val = Mpp.values;
  Real * Lpp_val = Lpp.values;

  UInt tangent_size = getTangentStiffnessVoigtSize<_bernoulli_beam_2>();
  UInt bt_d_b_size  = nb_nodes_per_element * nb_degre_of_freedom;
  b.clear();
  Vector<Real>::iterator<types::Matrix> B = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      (*B)(0,0) = Np_val[0];
      (*B)(0,1) = Np_val[1];
      (*B)(0,3) = Np_val[2];
      (*B)(0,4) = Np_val[3];

      (*B)(1,0) = Mpp_val[0];
      (*B)(1,1) = Mpp_val[1];
      (*B)(1,2) = Lpp_val[0];
      (*B)(1,3) = Mpp_val[2];
      (*B)(1,4) = Mpp_val[3];
      (*B)(1,5) = Lpp_val[1];
      ++B;

      Np_val  += 2*nb_nodes_per_element;
      Mpp_val += 2*nb_nodes_per_element;
      Lpp_val += 2*nb_nodes_per_element;
    }
  }
}

__END_AKANTU__
