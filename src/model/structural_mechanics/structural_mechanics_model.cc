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

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::StructuralMechanicsModel(Mesh & mesh,
						   UInt dim,
						   const ID & id,
						   const MemoryID & memory_id) :
  Model(mesh, dim, id, memory_id), Dumpable(),
  stress("stress", id, memory_id),
  element_material("element_material", id, memory_id),
  set_ID("beam sets", id, memory_id),
  stiffness_matrix(NULL),
  jacobian_matrix(NULL),
  solver(NULL),
  rotation_matrix("rotation_matices", id, memory_id) {
  AKANTU_DEBUG_IN();

  registerFEMObject<MyFEMType>("StructuralMechanicsFEM", mesh, spatial_dimension);

  this->displacement_rotation = NULL;
  this->force_momentum        = NULL;
  this->residual              = NULL;
  this->boundary              = NULL;
  this->increment             = NULL;

  if(spatial_dimension == 2)
    nb_degree_of_freedom = 3;
  else if (spatial_dimension == 3)
    nb_degree_of_freedom = 6;
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  this->registerDumper<DumperParaview>("paraview_all", id, true);
  this->addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_structural);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::~StructuralMechanicsModel() {
  AKANTU_DEBUG_IN();

  delete solver;
  delete stiffness_matrix;
  delete jacobian_matrix;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
/* -------------------------------------------------------------------------- */

void StructuralMechanicsModel::initFull(std::string material) {
  initModel();
  initArrays();
  initSolver();

  displacement_rotation->clear();
  force_momentum       ->clear();
  residual             ->clear();
  boundary             ->clear();
  increment            ->clear();


  Mesh::type_iterator it = getFEM().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEM().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    computeRotationMatrix(*it);
  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEM().getMesh().getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":boundary";
  std::stringstream sstr_incr; sstr_incr << id << ":increment";

  displacement_rotation = &(alloc<Real>(sstr_disp.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  force_momentum        = &(alloc<Real>(sstr_forc.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  residual              = &(alloc<Real>(sstr_resi.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  boundary              = &(alloc<bool>(sstr_boun.str(), nb_nodes, nb_degree_of_freedom, false));
  increment             = &(alloc<Real>(sstr_incr.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));

  Mesh::type_iterator it = getFEM().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEM().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    UInt nb_element = getFEM().getMesh().getNbElement(*it);
    UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(*it);

    element_material.alloc(nb_element, 1, *it, _not_ghost);
    set_ID.alloc(nb_element, 1, *it, _not_ghost);

    UInt size = getTangentStiffnessVoigtSize(*it);
    stress.alloc(nb_element * nb_quadrature_points, size , *it, _not_ghost);
  }

  dof_synchronizer = new DOFSynchronizer(getFEM().getMesh(), nb_degree_of_freedom);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initModel() {
  getFEM().initShapeFunctions(_not_ghost);
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initSolver(__attribute__((unused)) SolverOptions & options) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = getFEM().getMesh();
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_node = mesh.getNbGlobalNodes();

  std::stringstream sstr; sstr << id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(nb_global_node * nb_degree_of_freedom, _symmetric,
				     nb_degree_of_freedom, sstr.str(), memory_id);

  jacobian_matrix->buildProfile(mesh, *dof_synchronizer);

  std::stringstream sstr_sti; sstr_sti << id << ":stiffness_matrix";
  stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr_sti.str(), memory_id);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << id << ":solver";
  solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  solver->initialize(options);
#endif //AKANTU_HAS_SOLVER

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt StructuralMechanicsModel::getTangentStiffnessVoigtSize(const ElementType & type) {
  UInt size;
#define GET_TANGENT_STIFFNESS_VOIGT_SIZE(type)	\
  size = getTangentStiffnessVoigtSize<type>();

  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(GET_TANGENT_STIFFNESS_VOIGT_SIZE);
#undef GET_TANGENT_STIFFNESS_VOIGT_SIZE

  return size;
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  stiffness_matrix->clear();

  Mesh::type_iterator it = getFEM().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEM().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    ElementType type = *it;

#define ASSEMBLE_STIFFNESS_MATRIX(type)		\
    assembleStiffnessMatrix<type>();

    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(ASSEMBLE_STIFFNESS_MATRIX);
#undef ASSEMBLE_STIFFNESS_MATRIX
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeRotationMatrix<_bernoulli_beam_2>(Array<Real> & rotations){

  ElementType type = _bernoulli_beam_2;
  Mesh & mesh = getFEM().getMesh();
  UInt nb_element = mesh.getNbElement(type);

  Array<UInt>::iterator< Vector<UInt> > connec_it = mesh.getConnectivity(type).begin(2);
  Array<Real>::vector_iterator nodes_it = mesh.getNodes().begin(spatial_dimension);
  Array<Real>::matrix_iterator R_it = rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);

  for (UInt e = 0; e < nb_element; ++e, ++R_it, ++connec_it) {
    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x2;
    x2 = nodes_it[connec(1)]; // X2
    Vector<Real> x1;
    x1 = nodes_it[connec(0)]; // X1

    Real le = x1.distance(x2);
    Real c = (x2(0) - x1(0)) / le;
    Real s = (x2(1) - x1(1)) / le;

    /// Definition of the rotation matrix
    R(0,0) =  c;  R(0,1) = s;  R(0,2) = 0.;
    R(1,0) = -s;  R(1,1) = c;  R(1,2) = 0.;
    R(2,0) =  0.; R(2,1) = 0.; R(2,2) = 1.;
  }
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeRotationMatrix<_bernoulli_beam_3>(Array<Real> & rotations){
  ElementType type = _bernoulli_beam_3;
  Mesh & mesh = getFEM().getMesh();
  UInt nb_element = mesh.getNbElement(type);

  Array<Real>::vector_iterator n_it = mesh.getNormals(type).begin(spatial_dimension);
  Array<UInt>::iterator< Vector<UInt> > connec_it = mesh.getConnectivity(type).begin(2);
  Array<Real>::vector_iterator nodes_it = mesh.getNodes().begin(spatial_dimension);

  Matrix<Real> Pe    (spatial_dimension, spatial_dimension);
  Matrix<Real> Pg    (spatial_dimension, spatial_dimension);
  Matrix<Real> inv_Pg(spatial_dimension, spatial_dimension);
  Vector<Real> x_n(spatial_dimension); // x vect n

  Array<Real>::matrix_iterator R_it =
    rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);

  for (UInt e=0 ; e < nb_element; ++e, ++n_it, ++connec_it, ++R_it) {
    Vector<Real> & n = *n_it;
    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x;
    x = nodes_it[connec(1)]; // X2
    Vector<Real> y;
    y = nodes_it[connec(0)]; // X1

    Real l = x.distance(y);
    x -= y; // X2 - X1
    x_n.crossProduct(x, n);

    Pe.eye();
    Pe(0, 0) *=  l;
    Pe(1, 1) *= -l;

    Pg(0,0) = x(0); Pg(0,1) = x_n(0); Pg(0,2) = n(0);
    Pg(1,0) = x(1); Pg(1,1) = x_n(1); Pg(1,2) = n(1);
    Pg(2,0) = x(2); Pg(2,1) = x_n(2); Pg(2,2) = n(2);

    inv_Pg.inverse(Pg);

    Pe *= inv_Pg;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < spatial_dimension; ++j) {
	R(i, j) = Pe(i, j);
	R(i + spatial_dimension,j + spatial_dimension) = Pe(i, j);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeRotationMatrix(const ElementType & type) {
  Mesh & mesh = getFEM().getMesh();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_element = mesh.getNbElement(type);

  if(!rotation_matrix.exists(type)) {
    rotation_matrix.alloc(nb_element,
			  nb_degree_of_freedom*nb_nodes_per_element * nb_degree_of_freedom*nb_nodes_per_element,
			  type);
  } else {
    rotation_matrix(type).resize(nb_element);
  }
  rotation_matrix(type).clear();

  Array<Real>rotations(nb_element, nb_degree_of_freedom * nb_degree_of_freedom);
  rotations.clear();

#define COMPUTE_ROTATION_MATRIX(type)	\
  computeRotationMatrix<type>(rotations);

  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(COMPUTE_ROTATION_MATRIX);
#undef COMPUTE_ROTATION_MATRIX


  Array<Real>::matrix_iterator R_it = rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);
  Array<Real>::matrix_iterator T_it =
    rotation_matrix(type).begin(nb_degree_of_freedom*nb_nodes_per_element,
				nb_degree_of_freedom*nb_nodes_per_element);

  for (UInt el = 0; el < nb_element; ++el, ++R_it, ++T_it) {
    Matrix<Real> & T = *T_it;
    Matrix<Real> & R = *R_it;
    T.clear();
    for (UInt k = 0; k < nb_nodes_per_element; ++k){
      for (UInt i = 0; i < nb_degree_of_freedom; ++i)
	for (UInt j = 0; j < nb_degree_of_freedom; ++j)
	  T(k*nb_degree_of_freedom + i, k*nb_degree_of_freedom + j) = R(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeStresses() {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it = getFEM().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEM().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    ElementType type = *it;

#define COMPUTE_STRESS_ON_QUAD(type)		\
    computeStressOnQuad<type>();

    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(COMPUTE_STRESS_ON_QUAD);
#undef COMPUTE_STRESS_ON_QUAD
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::updateResidual() {
 AKANTU_DEBUG_IN();
 residual->copy(*force_momentum);

 Array<Real> ku(*displacement_rotation, true);
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
  jacobian_matrix->copyContent(*stiffness_matrix);
  jacobian_matrix->applyBoundary(*boundary);

  solver->setRHS(*residual);

  solver->solve(*increment);

  Real * increment_val     = increment->storage();
  Real * displacement_val  = displacement_rotation->storage();
  bool * boundary_val      = boundary->storage();

  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
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
  UInt nb_degree_of_freedom = displacement_rotation->getNbComponent();

  Real norm = 0;
  Real * increment_val     = increment->storage();
  bool * boundary_val      = boundary->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      if(!(*boundary_val) && is_local_node) {
	norm += *increment_val * *increment_val;
      }
      boundary_val++;
      increment_val++;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&norm, 1, _so_sum);

  error = sqrt(norm);
  AKANTU_DEBUG_ASSERT(!isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (error < tolerance);
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_2>(Array<Real> & tangent_moduli) {
  UInt nb_element                 = getFEM().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(_bernoulli_beam_2);
  UInt tangent_size = 2;

  Array<Real>::matrix_iterator D_it = tangent_moduli.begin(tangent_size, tangent_size);
  Array<UInt> & el_mat = element_material(_bernoulli_beam_2, _not_ghost);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = el_mat(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real I = materials[mat].I;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0,0) = E * A;
      D(1,1) = E * I;
    }
  }
}
/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_3>(Array<Real> & tangent_moduli) {
  UInt nb_element                 = getFEM().getMesh().getNbElement(_bernoulli_beam_3);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(_bernoulli_beam_3);
  UInt tangent_size = 4;

  Array<Real>::matrix_iterator D_it = tangent_moduli.begin(tangent_size, tangent_size);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = element_material(_bernoulli_beam_3, _not_ghost)(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real Iz = materials[mat].Iz;
    Real Iy = materials[mat].Iy;
    Real GJ = materials[mat].GJ;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0,0) = E * A;
      D(1,1) = E * Iz;
      D(2,2) = E * Iy;
      D(3,3) = GJ;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<_bernoulli_beam_2>(Array<Real> & b, bool local) {
  UInt nb_element                 = getFEM().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(_bernoulli_beam_2);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(_bernoulli_beam_2);

  MyFEMType & fem = getFEMClass<MyFEMType>();
  Array<Real>::const_vector_iterator shape_Np  = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Mpp = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Lpp = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 2).begin(nb_nodes_per_element);

  UInt tangent_size = getTangentStiffnessVoigtSize<_bernoulli_beam_2>();
  UInt bt_d_b_size  = nb_nodes_per_element * nb_degree_of_freedom;
  b.clear();
  Array<Real>::matrix_iterator B_it = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & B = *B_it;
      const Vector<Real> & Np  = *shape_Np;
      const Vector<Real> & Lpp = *shape_Lpp;
      const Vector<Real> & Mpp = *shape_Mpp;

      B(0,0) = Np(0);
      B(0,3) = Np(1);

      B(1,1) = Mpp(0);
      B(1,2) = Lpp(0);
      B(1,4) = Mpp(1);
      B(1,5) = Lpp(1);

      ++B_it;
      ++shape_Np;
      ++shape_Mpp;
      ++shape_Lpp;
    }

    // ++R_it;
  }
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<_bernoulli_beam_3>(Array<Real> & b, bool local) {
  MyFEMType & fem = getFEMClass<MyFEMType>();

  UInt nb_element                 = getFEM().getMesh().getNbElement(_bernoulli_beam_3);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(_bernoulli_beam_3);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(_bernoulli_beam_3);

  Array<Real>::const_vector_iterator shape_Np  = fem.getShapesDerivatives(_bernoulli_beam_3, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Mpp = fem.getShapesDerivatives(_bernoulli_beam_3, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Lpp = fem.getShapesDerivatives(_bernoulli_beam_3, _not_ghost, 2).begin(nb_nodes_per_element);

  UInt tangent_size = getTangentStiffnessVoigtSize<_bernoulli_beam_3>();
  UInt bt_d_b_size  = nb_nodes_per_element * nb_degree_of_freedom;

  b.clear();

  Array<Real>::matrix_iterator B_it = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & B  = *B_it;

      const Vector<Real> & Np  = *shape_Np;
      const Vector<Real> & Lpp = *shape_Lpp;
      const Vector<Real> & Mpp = *shape_Mpp;

      B(0,0)  =  Np(0);
      B(0,6)  =  Np(1);

      B(1,1)  =  Mpp(0);
      B(1,5)  =  Lpp(0);
      B(1,7)  =  Mpp(1);
      B(1,11) =  Lpp(1);

      B(2,2)  =  Mpp(0);
      B(2,4)  = -Lpp(0);
      B(2,8)  =  Mpp(1);
      B(2,10) = -Lpp(1);

      B(3,3)  =  Np(0);
      B(3,9)  =  Np(1);

      ++B_it;
      ++shape_Np;
      ++shape_Mpp;
      ++shape_Lpp;
    }
  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::addDumpFieldToDumper(const std::string & dumper_name,
						    const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
#define ADD_FIELD(dumper_name, id, field, type, n, stride)		\
  internalAddDumpFieldToDumper(dumper_name, id,				\
		       new DumperIOHelper::NodalField<type>(*field, n, stride))

  UInt n;
  if(spatial_dimension == 2) {
    n = 2;
  } else n = 3;

  if(field_id == "displacement" ) { ADD_FIELD(dumper_name, "displacement", displacement_rotation, Real,
						  n, 0); }
  else if(field_id == "rotation") { ADD_FIELD(dumper_name, "rotation", displacement_rotation, Real,
						  nb_degree_of_freedom - n, n); }
  else if(field_id == "force"   ) { ADD_FIELD(dumper_name, "force", force_momentum, Real,
						  n, 0); }
  else if(field_id == "momentum") { ADD_FIELD(dumper_name, "momentum", force_momentum, Real,
						  nb_degree_of_freedom - n, n); }
  else if(field_id == "residual") { ADD_FIELD(dumper_name, "residual", residual, Real, nb_degree_of_freedom, 0); }
  else if(field_id == "boundary") { ADD_FIELD(dumper_name, "boundary", boundary, bool, nb_degree_of_freedom, 0); }
  else if(field_id == "element_index_by_material") {
    internalAddDumpFieldToDumper(dumper_name, 
				 field_id,
			 new DumperIOHelper::ElementalField<UInt>(element_material,
								  spatial_dimension,
								  _not_ghost,
								  _ek_regular));
  }
#undef ADD_FIELD
#endif
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::addDumpFieldVectorToDumper(const std::string & dumper_name,
							  const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
#define ADD_FIELD(dumper_name, id, field, type, n, stride)		\
  DumperIOHelper::Field * f =						\
    new DumperIOHelper::NodalField<type>(*field, n, stride);	\
  f->setPadding(3);							\
  internalAddDumpFieldToDumper(dumper_name, id, f)

  UInt n;
  if(spatial_dimension == 2) {
    n = 2;
  } else n = 3;

  if(field_id == "displacement" ) { ADD_FIELD(dumper_name, "displacement", displacement_rotation, Real,
					      n, 0); }
  else if(field_id == "force"   ) { ADD_FIELD(dumper_name, "force", force_momentum, Real,
					      n, 0); }
#undef ADD_FIELD
#endif
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::addDumpFieldTensorToDumper(__attribute__((unused)) const std::string & dumper_name,
							  __attribute__((unused)) const std::string & field_id) {
}


__END_AKANTU__
