/**
 * @file   solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 14:35:38 2010
 *
 * @brief  Implementation of the SolidMechanicsModel class
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
#include "solid_mechanics_model.hh"
#include "aka_math.hh"
#include "integration_scheme_2nd_order.hh"

#include "static_communicator.hh"

#include "dof_synchronizer.hh"


#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
/**
 * A solid mechanics model need a mesh  and a dimension to be created. the model
 * by it  self can not  do a lot,  the good init  functions should be  called in
 * order to configure the model depending on what we want to do.
 *
 * @param  mesh mesh  representing  the model  we  want to  simulate
 * @param dim spatial  dimension of the problem, if dim =  0 (default value) the
 * dimension of the problem is assumed to be the on of the mesh
 * @param id an id to identify the model
 */
SolidMechanicsModel::SolidMechanicsModel(Mesh & mesh,
					 UInt dim,
					 const ID & id,
					 const MemoryID & memory_id) :
  Model(id, memory_id),
  time_step(NAN), f_m2a(1.0),
  mass_matrix(NULL),
  velocity_damping_matrix(NULL),
  stiffness_matrix(NULL),
  jacobian_matrix(NULL),
  element_material("element_material", id),
  integrator(NULL),
  increment_flag(false), solver(NULL),
  spatial_dimension(dim), mesh(mesh), dynamic(false), implicit(false) {
  AKANTU_DEBUG_IN();

  createSynchronizerRegistry(this);

  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();
  registerFEMObject<MyFEMType>("SolidMechanicsFEM", mesh, spatial_dimension);

  this->displacement = NULL;
  this->mass         = NULL;
  this->velocity     = NULL;
  this->acceleration = NULL;
  this->force        = NULL;
  this->residual     = NULL;
  this->boundary     = NULL;

  this->increment    = NULL;
  this->increment_acceleration = NULL;

  this->dof_synchronizer = NULL;

  materials.clear();

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
SolidMechanicsModel::~SolidMechanicsModel() {
  AKANTU_DEBUG_IN();

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    delete *mat_it;
  }
  materials.clear();

  delete integrator;

  if(solver) delete solver;
  if(mass_matrix) delete mass_matrix;
  if(velocity_damping_matrix) delete velocity_damping_matrix;

  if(stiffness_matrix && stiffness_matrix != jacobian_matrix)
    delete stiffness_matrix;

  if(jacobian_matrix) delete jacobian_matrix;

  if(dof_synchronizer) delete dof_synchronizer;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
/* -------------------------------------------------------------------------- */
/**
 * This function groups  many of the initialization in on  function. For most of
 * basics  case the  function should  be  enough. The  functions initialize  the
 * model, the internal  vectors, set them to 0, and  depending on the parameters
 * it also initialize the explicit or implicit solver.
 *
 * @param material_file the file containing the materials to use
 * @param implicit_scheme set to true for implicit simulations, or to false for
 * explicit ones
 * @param implicit_dynamic in case of implicit solver activated it specify if we
 * want to solve a static or dynamic case.
 */
void SolidMechanicsModel::initFull(std::string material_file,
				   bool implicit_scheme,
				   bool implicit_dynamic) {
  // initialize the model
  initModel();

  // initialize the vectors
  initVectors();

  // set the initial condition to 0
  force->clear();
  velocity->clear();
  acceleration->clear();
  displacement->clear();

  // initialize pcb
  if(pbc_pair.size()!=0)
    initPBC();

  // initialize the time integration schemes
  if(implicit_scheme)
    initImplicit(implicit_dynamic);
  else
    initExplicit();

  // initialize the materials
  if(material_file != "") {
    readMaterials(material_file);
    initMaterials();
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initParallel(MeshPartition * partition,
				       DataAccessor * data_accessor) {
  AKANTU_DEBUG_IN();

  if (data_accessor == NULL) data_accessor = this;
  Synchronizer & synch_parallel = createParallelSynch(partition,data_accessor);

  synch_registry->registerSynchronizer(synch_parallel,_gst_smm_mass);
  synch_registry->registerSynchronizer(synch_parallel,_gst_smm_for_strain);
  synch_registry->registerSynchronizer(synch_parallel,_gst_smm_boundary);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initCohesive() {
  registerFEMObject<MyFEMCohesiveType>("CohesiveFEM", mesh, spatial_dimension);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initFEMBoundary(bool create_surface) {

  if(create_surface)
    MeshUtils::buildFacets(mesh, true, false);

  FEM & fem_boundary = getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initExplicit() {
  AKANTU_DEBUG_IN();

  if (integrator) delete integrator;
  integrator = new CentralDifference();

  increment_acceleration = new Vector<Real>(this->mesh.getNbNodes(), 
					    this->spatial_dimension, 0.);

  dynamic = true;
  implicit = false;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/**
 * Allocate all the needed vectors. By  default their are not necessarily set to
 * 0
 *
 */
void SolidMechanicsModel::initVectors() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  std::stringstream sstr_mass; sstr_mass << id << ":mass";
  std::stringstream sstr_velo; sstr_velo << id << ":velocity";
  std::stringstream sstr_acce; sstr_acce << id << ":acceleration";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":boundary";

  displacement = &(alloc<Real>(sstr_disp.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  mass         = &(alloc<Real>(sstr_mass.str(), nb_nodes, spatial_dimension, 0));
  velocity     = &(alloc<Real>(sstr_velo.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  acceleration = &(alloc<Real>(sstr_acce.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  force        = &(alloc<Real>(sstr_forc.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  residual     = &(alloc<Real>(sstr_resi.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  boundary     = &(alloc<bool>(sstr_boun.str(), nb_nodes, spatial_dimension, false));

  std::stringstream sstr_curp; sstr_curp << id << ":current_position";
  current_position = &(alloc<Real>(sstr_curp.str(), 0, spatial_dimension, REAL_INIT_VALUE));

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      element_material.alloc(nb_element, 1, *it, gt);
    }
  }

  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 *
 */
void SolidMechanicsModel::initModel() {
  /// \todo add  the current position  as a parameter to  initShapeFunctions for
  /// large deformation
  getFEM().initShapeFunctions(_not_ghost);
  getFEM().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initPBC() {
  Model::initPBC();
  registerPBCSynchronizer();

  // as long as there are ones on the diagonal of the matrix, we can put boudandary true for slaves
  std::map<UInt, UInt>::iterator it = pbc_pair.begin();
  std::map<UInt, UInt>::iterator end = pbc_pair.end();
  UInt dim = mesh.getSpatialDimension();
  while(it != end) {
    for (UInt i=0; i<dim; ++i)
      (*boundary)((*it).first,i) = true;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::registerPBCSynchronizer(){
  PBCSynchronizer * synch = new PBCSynchronizer(pbc_pair);
  synch_registry->registerSynchronizer(*synch, _gst_smm_uv);
  synch_registry->registerSynchronizer(*synch, _gst_smm_mass);
  changeLocalEquationNumberforPBC(pbc_pair, mesh.getSpatialDimension());
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  current_position->resize(nb_nodes);
  //Vector<Real> * current_position = new Vector<Real>(nb_nodes, spatial_dimension, NAN, "position");
  Real * current_position_val = current_position->values;
  Real * position_val         = mesh.getNodes().values;
  Real * displacement_val     = displacement->values;

  /// compute current_position = initial_position + displacement
  memcpy(current_position_val, position_val, nb_nodes*spatial_dimension*sizeof(Real));
  for (UInt n = 0; n < nb_nodes*spatial_dimension; ++n) {
    *current_position_val++ += *displacement_val++;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initializeUpdateResidualData() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();
  residual->resize(nb_nodes);

  /// copy the forces in residual for boundary conditions
  memcpy(residual->values, force->values, nb_nodes*spatial_dimension*sizeof(Real));

  updateCurrentPosition();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* Explicit scheme                                                            */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/**
 * This function compute  the second member of the motion  equation.  That is to
 * say the sum of forces @f$ r = F_{ext} - F_{int} @f$. @f$ F_{ext} @f$ is given
 * by the  user in  the force  vector, and @f$  F_{int} @f$  is computed  as @f$
 * F_{int} = \int_{\Omega} N \sigma d\Omega@f$
 *
 */
void SolidMechanicsModel::updateResidual(bool need_initialize) {
  AKANTU_DEBUG_IN();

  // f = f_ext - f_int

  // f = f_ext
  if (need_initialize) initializeUpdateResidualData();

  // f -= fint
  // start synchronization
  synch_registry->asynchronousSynchronize(_gst_smm_uv);
  synch_registry->waitEndSynchronize(_gst_smm_uv);
  synch_registry->asynchronousSynchronize(_gst_smm_for_strain);

  // call update residual on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    //(*mat_it)->updateResidual(*current_position, _not_ghost);
    (*mat_it)->updateResidual(*displacement, _not_ghost);
  }

  // finalize communications
  synch_registry->waitEndSynchronize(_gst_smm_for_strain);

  // call update residual on each ghost elements
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    //(*mat_it)->updateResidual(*current_position, _ghost);
    (*mat_it)->updateResidual(*displacement, _ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateResidualInternal() {
  AKANTU_DEBUG_IN();

  // f = f_ext - f_int - Ma - Cv = r - Ma - Cv;

  if(dynamic) {
    // f -= Ma
    if(mass_matrix) {
      // if full mass_matrix
      Vector<Real> * Ma = new Vector<Real>(*acceleration, true, "Ma");
      *Ma *= *mass_matrix;
      *residual -= *Ma;
      delete Ma;
    } else {
      // else lumped mass
      UInt nb_nodes = acceleration->getSize();
      UInt nb_degree_of_freedom = acceleration->getNbComponent();

      Real * mass_val     = mass->values;
      Real * accel_val    = acceleration->values;
      Real * res_val      = residual->values;
      bool * boundary_val = boundary->values;

      for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
	if(!(*boundary_val)) {
	  *res_val -= *accel_val * *mass_val;
	}
	boundary_val++;
	res_val++;
	mass_val++;
	accel_val++;
      }
    }

    // f -= Cv
    if(velocity_damping_matrix) {
      Vector<Real> * Cv = new Vector<Real>(*velocity);
      *Cv *= *velocity_damping_matrix;
      *residual -= *Cv;
      delete Cv;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateAcceleration() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = acceleration->getSize();
  UInt nb_degree_of_freedom = acceleration->getNbComponent();

  updateResidualInternal();

  if(!implicit && !mass_matrix) {

    Real * mass_val     = mass->values;
    Real * residual_val = residual->values;
    bool * boundary_val = boundary->values;
    Real * inc = increment_acceleration->values;
    Real * accel_val    = acceleration->values;

    for (UInt n = 0; n < nb_nodes; ++n) {
      for (UInt d = 0; d < nb_degree_of_freedom; d++) {
	if(!(*boundary_val)) {
	  *inc = f_m2a * (*residual_val / *mass_val);
	}
	residual_val++;
	boundary_val++;
	inc++;
	mass_val++;
	accel_val++;
      }
    }
  } else {
    solveDynamic<NewmarkBeta::_acceleration_corrector>(*increment_acceleration);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitPred() {
  AKANTU_DEBUG_IN();

  if(increment_flag) {
    memcpy(increment->values,
	   displacement->values,
	   displacement->getSize()*displacement->getNbComponent()*sizeof(Real));
  }

  AKANTU_DEBUG_ASSERT(integrator,"itegrator should have been allocated: "
		      << "have called initExplicit ? "
		      << "or initImplicit ?");

  integrator->integrationSchemePred(time_step,
				    *displacement,
				    *velocity,
				    *acceleration,
				    *boundary);

  if(increment_flag) {
    Real * inc_val = increment->values;
    Real * dis_val = displacement->values;
    UInt nb_nodes = displacement->getSize();

    for (UInt n = 0; n < nb_nodes; ++n) {
      *inc_val = *dis_val - *inc_val;
      inc_val++;
      dis_val++;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitCorr() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemeCorrAccel(time_step,
					 *displacement,
					 *velocity,
					 *acceleration,
					 *boundary,
					 *increment_acceleration);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Implicit scheme                                                            */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/**
 * Initialize the solver and create the sparse matrices needed.
 *
 */
void SolidMechanicsModel::initSolver(SolverOptions & options) {
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_node = mesh.getNbGlobalNodes();

  std::stringstream sstr; sstr << id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(nb_global_node * spatial_dimension, _symmetric,
					 spatial_dimension, sstr.str(), memory_id);
  //  dof_synchronizer->initGlobalDOFEquationNumbers();
  jacobian_matrix->buildProfile(mesh, *dof_synchronizer);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << id << ":solver";
  solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  solver->initialize(options);
#endif //AKANTU_HAS_SOLVER
}


/* -------------------------------------------------------------------------- */
/**
 * Initialize the implicit solver, either for dynamic or static cases, 
 *
 * @param dynamic
 */
void SolidMechanicsModel::initImplicit(bool dynamic, SolverOptions & solver_options) {
  AKANTU_DEBUG_IN();

  this->dynamic = dynamic;
  implicit = true;

  initSolver(solver_options);

  if(dynamic) {
    if(integrator) delete integrator;
    integrator = new TrapezoidalRule2();

    std::stringstream sstr; sstr << id << ":stiffness_matrix";
    stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr.str(), memory_id);
  } else {
    stiffness_matrix = jacobian_matrix;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initialAcceleration() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving Ma = f");

  Solver * acc_solver = NULL;
#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr; sstr << id << ":solver_mass_matrix";
  acc_solver = new SolverMumps(*mass_matrix, sstr.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS
  acc_solver->initialize();

  mass_matrix->applyBoundary(*boundary);

  acc_solver->setRHS(*residual);
  acc_solver->solve(*acceleration);

  delete acc_solver;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  stiffness_matrix->clear();

  // call compute stiffness matrix on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->assembleStiffnessMatrix(*displacement, _not_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::solveDynamic() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
		      "You should first initialize the implicit solver and assemble the stiffness matrix");
  AKANTU_DEBUG_ASSERT(mass_matrix != NULL,
		      "You should first initialize the implicit solver and assemble the mass matrix");

  if(!increment) setIncrementFlagOn();

  updateResidualInternal();

  solveDynamic<NewmarkBeta::_displacement_corrector>(*increment);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::solveStatic() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving Ku = f");
  AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
		      "You should first initialize the implicit solver and assemble the stiffness matrix");

  UInt nb_nodes = displacement->getSize();
  UInt nb_degree_of_freedom = displacement->getNbComponent();

  stiffness_matrix->applyBoundary(*boundary);

  if(dynamic)
    jacobian_matrix->copyContent(*stiffness_matrix);

  solver->setRHS(*residual);

  if(!increment) setIncrementFlagOn();

  solver->solve(*increment);

  Real * increment_val     = increment->values;
  Real * displacement_val  = displacement->values;
  bool * boundary_val      = boundary->values;

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
bool SolidMechanicsModel::testConvergenceIncrement(Real tolerance) {
  Real error;
  bool tmp = testConvergenceIncrement(tolerance, error);

  AKANTU_DEBUG_INFO("Norm of increment : " << error);

  return tmp;
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceIncrement(Real tolerance, Real & error) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = displacement->getSize();
  UInt nb_degree_of_freedom = displacement->getNbComponent();

  error = 0;
  Real norm[2] = {0., 0.};
  Real * increment_val    = increment->storage();
  bool * boundary_val     = boundary->storage();
  Real * displacement_val = displacement->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      if(!(*boundary_val) && is_local_node) {
	norm[0] += *increment_val * *increment_val;
	norm[1] += *displacement_val * *displacement_val;
      }
      boundary_val++;
      increment_val++;
      displacement_val++;
    }
  }

  StaticCommunicator::getStaticCommunicator()->allReduce(norm, 2, _so_sum);

  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);
  AKANTU_DEBUG_ASSERT(!Math::isnan(norm[0]), "Something goes wrong in the solve phase");


  AKANTU_DEBUG_OUT();
  error = norm[0] / norm[1];
  return (error < tolerance);
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceResidual(Real tolerance) {
  Real error;
  bool tmp = testConvergenceResidual(tolerance, error);

  AKANTU_DEBUG_INFO("Norm of residual : " << error);

  return tmp;
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceResidual(Real tolerance, Real & norm) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = residual->getSize();

  norm = 0;
  Real * residual_val = residual->values;
  bool * boundary_val = boundary->values;

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    if(is_local_node) {
      for (UInt d = 0; d < spatial_dimension; ++d) {
	if(!(*boundary_val)) {
	  norm += *residual_val * *residual_val;
	}
	boundary_val++;
	residual_val++;
      }
    } else {
      boundary_val += spatial_dimension;
      residual_val += spatial_dimension;
    }
  }

  StaticCommunicator::getStaticCommunicator()->allReduce(&norm, 1, _so_sum);

  norm = sqrt(norm);

  AKANTU_DEBUG_ASSERT(!Math::isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (norm < tolerance);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::implicitPred() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemePred(time_step,
				    *displacement,
				    *velocity,
				    *acceleration,
				    *boundary);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::implicitCorr() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemeCorrDispl(time_step,
					 *displacement,
					 *velocity,
					 *acceleration,
					 *boundary,
					 *increment);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* Information                                                                */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::synchronizeBoundaries() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(synch_registry,"Synchronizer registry was not initialized."
		      << " Did you call initParallel ?");
  synch_registry->synchronize(_gst_smm_boundary);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setIncrementFlagOn() {
  AKANTU_DEBUG_IN();

  if(!increment) {
    UInt nb_nodes = mesh.getNbNodes();
    std::stringstream sstr_inc; sstr_inc << id << ":increment";
    increment = &(alloc<Real>(sstr_inc.str(), nb_nodes, spatial_dimension, 0));
  }

  increment_flag = true;

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real min_dt = getStableTimeStep(_not_ghost);

  /// reduction min over all processors
  StaticCommunicator::getStaticCommunicator()->allReduce(&min_dt, 1, _so_min);

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));
  Real min_dt = std::numeric_limits<Real>::max();

  Real * coord    = mesh.getNodes().values;
  Real * disp_val = displacement->values;

  Element elem;
  elem.ghost_type = ghost_type;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != end; ++it) {
    elem.type = *it;
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
    UInt nb_element           = mesh.getNbElement(*it);

    UInt * conn         = mesh.getConnectivity(*it, ghost_type).storage();
    UInt * elem_mat_val = element_material(*it, ghost_type).storage();

    Real * u = new Real[nb_nodes_per_element*spatial_dimension];

    for (UInt el = 0; el < nb_element; ++el) {
      UInt el_offset  = el * nb_nodes_per_element;
      elem.element = el;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	UInt offset_conn = conn[el_offset + n] * spatial_dimension;
	memcpy(u + n * spatial_dimension,
	       coord + offset_conn,
	       spatial_dimension * sizeof(Real));

	for (UInt i = 0; i < spatial_dimension; ++i) {
	  u[n * spatial_dimension + i] += disp_val[offset_conn + i];
	}
      }

      Real el_size    = getFEM().getElementInradius(u, *it);
      Real el_dt      = mat_val[elem_mat_val[el]]->getStableTimeStep(el_size, elem);
      min_dt = min_dt > el_dt ? el_dt : min_dt;
    }

    delete [] u;
  }

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  /// call update residual on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    epot += (*mat_it)->getPotentialEnergy();
  }

  /// reduction sum over all processors
  StaticCommunicator::getStaticCommunicator()->allReduce(&epot, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy() {
  AKANTU_DEBUG_IN();

  Real ekin = 0.;

  UInt nb_nodes = mesh.getNbNodes();

  Real * vel_val  = velocity->values;
  Real * mass_val = mass->values;

  for (UInt n = 0; n < nb_nodes; ++n) {
    Real mv2 = 0;
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool is_not_pbc_slave_node = !getIsPBCSlaveNode(n);
    for (UInt i = 0; i < spatial_dimension; ++i) {
      //      if(is_local_node) {
      mv2 += is_local_node * is_not_pbc_slave_node * *vel_val * *vel_val * *mass_val;
      //      }
      vel_val++;
      mass_val++;
    }
    ekin += mv2;
  }

  StaticCommunicator::getStaticCommunicator()->allReduce(&ekin, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return ekin * .5;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Solid Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << spatial_dimension << std::endl;

  stream << space << " + fem [" << std::endl;
  getFEM().printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;
  stream << space << " + nodals information [" << std::endl;
  displacement->printself(stream, indent + 2);
  mass        ->printself(stream, indent + 2);
  velocity    ->printself(stream, indent + 2);
  acceleration->printself(stream, indent + 2);
  force       ->printself(stream, indent + 2);
  residual    ->printself(stream, indent + 2);
  boundary    ->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + connectivity type information [" << std::endl;
  element_material.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */


__END_AKANTU__
