/**
 * @file   solid_mechanics_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 27 18:15:37 2010
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
#include "aka_math.hh"
#include "aka_common.hh"
#include "solid_mechanics_model.hh"
#include "integration_scheme_2nd_order.hh"

#include "static_communicator.hh"

#include "dof_synchronizer.hh"
#include <cmath>

#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_iohelper_tmpl_homogenizing_field.hh"
#  include "dumper_iohelper_tmpl_material_internal_field.hh"
#endif

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
  Model(mesh, id, memory_id), Dumpable<DumperParaview>(id),
  time_step(NAN), f_m2a(1.0),
  mass_matrix(NULL),
  velocity_damping_matrix(NULL),
  stiffness_matrix(NULL),
  jacobian_matrix(NULL),
  //element_material("element_material", id),
  element_index_by_material("element index by material", id),
  integrator(NULL),
  increment_flag(false), solver(NULL),
  spatial_dimension(dim),
  synch_parallel(NULL) {
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

  mesh.registerEventHandler(*this);

  addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);

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

  if(synch_parallel) delete synch_parallel;

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
 * @param material_file the  file containing the materials to  use
 * @param method the analysis method wanted.  See the akantu::AnalysisMethod for
 * the different possibilites
 */
void SolidMechanicsModel::initFull(std::string material_file,
				   AnalysisMethod analysis_method) {
  method = analysis_method;

  // initialize the model
  initModel();

  // initialize the vectors
  initArrays();

  // set the initial condition to 0
  force->clear();
  velocity->clear();
  acceleration->clear();
  displacement->clear();

  // initialize pcb
  if(pbc_pair.size()!=0)
    initPBC();

  // initialize the time integration schemes
  switch(method) {
  case _explicit_dynamic:
    initExplicit();
    break;
  case _implicit_dynamic:
    initImplicit(true);
    break;
  case _static:
    initImplicit(false);
    break;
  }

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
  synch_parallel = &createParallelSynch(partition,data_accessor);

  synch_registry->registerSynchronizer(*synch_parallel, _gst_material_id);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_mass);
  //  synch_registry->registerSynchronizer(synch_parallel, _gst_smm_for_strain);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_stress);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_boundary);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initFEMBoundary(bool create_surface) {

  if(create_surface)
    MeshUtils::buildFacets(mesh);

  FEM & fem_boundary = getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initExplicit() {
  AKANTU_DEBUG_IN();

  method = _explicit_dynamic;

  if (integrator) delete integrator;
  integrator = new CentralDifference();

  UInt nb_nodes = acceleration->getSize();
  UInt nb_degree_of_freedom = acceleration->getNbComponent();

  std::stringstream sstr; sstr << id << ":increment_acceleration";
  increment_acceleration = &(alloc<Real>(sstr.str(), nb_nodes, nb_degree_of_freedom, Real()));

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/**
 * Allocate all the needed vectors. By  default their are not necessarily set to
 * 0
 *
 */
void SolidMechanicsModel::initArrays() {
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
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      element_index_by_material.alloc(nb_element, 2, *it, gt);
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
  synch_registry->registerSynchronizer(*synch, _gst_smm_res);
  changeLocalEquationNumberforPBC(pbc_pair, mesh.getSpatialDimension());
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  current_position->resize(nb_nodes);
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

  // start synchronization
  synch_registry->asynchronousSynchronize(_gst_smm_uv);
  synch_registry->waitEndSynchronize(_gst_smm_uv);

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

  // communicate the displacement
  // synch_registry->asynchronousSynchronize(_gst_smm_for_strain);

  std::vector<Material *>::iterator mat_it;

  // call update residual on each local elements
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.computeAllStresses(_not_ghost);
  }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  /* ------------------------------------------------------------------------ */
  /* Computation of the non local part */
  synch_registry->asynchronousSynchronize(_gst_mnl_for_average);

  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.computeAllNonLocalStresses(_not_ghost);
  }

  synch_registry->waitEndSynchronize(_gst_mnl_for_average);

  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.computeAllNonLocalStresses(_ghost);
  }
#endif

  /* ------------------------------------------------------------------------ */
  /* assembling the forces internal */
  // communicate the strain
  synch_registry->asynchronousSynchronize(_gst_smm_stress);

  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.assembleResidual(_not_ghost);
  }

  // finalize communications
  synch_registry->waitEndSynchronize(_gst_smm_stress);
  //  synch_registry->waitEndSynchronize(_gst_smm_for_strain);

  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.assembleResidual(_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeStresses() {
  if (method == _explicit_dynamic) {
    // start synchronization
    synch_registry->asynchronousSynchronize(_gst_smm_uv);
    synch_registry->waitEndSynchronize(_gst_smm_uv);

    // compute stresses on all local elements for each materials
    std::vector<Material *>::iterator mat_it;
    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllStresses(_not_ghost);
    }


    /* ------------------------------------------------------------------------ */
#ifdef AKANTU_DAMAGE_NON_LOCAL
    /* Computation of the non local part */
    synch_registry->asynchronousSynchronize(_gst_mnl_for_average);

    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllNonLocalStresses(_not_ghost);
    }

    synch_registry->waitEndSynchronize(_gst_mnl_for_average);

    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllNonLocalStresses(_ghost);
    }
#endif
  } else {
    std::vector<Material *>::iterator mat_it;
    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllStressesFromTangentModuli(_not_ghost);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateResidualInternal() {
  AKANTU_DEBUG_IN();

  // f = f_ext - f_int - Ma - Cv = r - Ma - Cv;

  if(method != _static) {
    // f -= Ma
    if(mass_matrix) {
      // if full mass_matrix
      Array<Real> * Ma = new Array<Real>(*acceleration, true, "Ma");
      *Ma *= *mass_matrix;
      /// \todo check unit conversion for implicit dynamics
      //      *Ma /= f_m2a
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
	  *res_val -= *accel_val * *mass_val /f_m2a;
	}
	boundary_val++;
	res_val++;
	mass_val++;
	accel_val++;
      }
    }

    // f -= Cv
    if(velocity_damping_matrix) {
      Array<Real> * Cv = new Array<Real>(*velocity);
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

  updateResidualInternal();

  if(method == _explicit_dynamic && !mass_matrix) {
    /* residual = residual_{n+1} - M * acceleration_n therefore
       solution = increment acceleration not acceleration */
    solveLumped(*increment_acceleration,
		*mass,
		*residual,
		*boundary,
		f_m2a);
  } else {
    solveDynamic<NewmarkBeta::_acceleration_corrector>(*increment_acceleration);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::solveLumped(Array<Real> & x,
				      const Array<Real> & A,
				      const Array<Real> & b,
				      const Array<bool> & boundary,
				      Real alpha) {
  Real * A_val = A.storage();
  Real * b_val = b.storage();
  Real * x_val = x.storage();
  bool * boundary_val = boundary.storage();

  UInt nb_degrees_of_freedom = x.getSize() * x.getNbComponent();

  for (UInt n = 0; n < nb_degrees_of_freedom; ++n) {
    if(!(*boundary_val)) {
      *x_val = alpha * (*b_val / *A_val);
    }
    x_val++;
    A_val++;
    b_val++;
    boundary_val++;
  }
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
void SolidMechanicsModel::initSolver(__attribute__((unused)) SolverOptions & options) {
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_node = mesh.getNbGlobalNodes();

  std::stringstream sstr; sstr << id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(nb_global_node * spatial_dimension, _symmetric,
					 spatial_dimension, sstr.str(), memory_id);

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
}


/* -------------------------------------------------------------------------- */
/**
 * Initialize the implicit solver, either for dynamic or static cases,
 *
 * @param dynamic
 */
void SolidMechanicsModel::initImplicit(bool dynamic, SolverOptions & solver_options) {
  AKANTU_DEBUG_IN();

  method = dynamic ? _implicit_dynamic : _static;

  initSolver(solver_options);

  if(method == _implicit_dynamic) {
    if(integrator) delete integrator;
    integrator = new TrapezoidalRule2();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initialAcceleration() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving Ma = f");

  Solver * acc_solver = NULL;

  std::stringstream sstr; sstr << id << ":tmp_mass_matrix";
  SparseMatrix * tmp_mass = new SparseMatrix(*mass_matrix, sstr.str(), memory_id);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solver; sstr << id << ":solver_mass_matrix";
  acc_solver = new SolverMumps(*mass_matrix, sstr_solver.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  acc_solver->initialize();

  tmp_mass->applyBoundary(*boundary);

  acc_solver->setRHS(*residual);
  acc_solver->solve(*acceleration);

  delete acc_solver;
  delete tmp_mass;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  stiffness_matrix->clear();

  // call compute stiffness matrix on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->assembleStiffnessMatrix(_not_ghost);
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

  //  if(method != _static)
  jacobian_matrix->copyContent(*stiffness_matrix);
  jacobian_matrix->applyBoundary(*boundary);

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

  StaticCommunicator::getStaticCommunicator().allReduce(norm, 2, _so_sum);

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

  StaticCommunicator::getStaticCommunicator().allReduce(&norm, 1, _so_sum);

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
		      << " Did you call initParallel?");
  synch_registry->synchronize(_gst_smm_boundary);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::synchronizeResidual() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(synch_registry,"Synchronizer registry was not initialized."
		      << " Did you call initPBC?");
  synch_registry->synchronize(_gst_smm_res);
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
  StaticCommunicator::getStaticCommunicator().allReduce(&min_dt, 1, _so_min);

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));
  Real min_dt = std::numeric_limits<Real>::max();

  updateCurrentPosition();

  Element elem;
  elem.ghost_type = ghost_type;
  elem.kind = _ek_regular;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != end; ++it) {
    elem.type = *it;
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
    UInt nb_element           = mesh.getNbElement(*it);

    Array<UInt>::iterator< Vector<UInt> > eibm =
      element_index_by_material(*it, ghost_type).begin(2);

    Array<Real> X(0, nb_nodes_per_element*spatial_dimension);
    FEM::extractNodalToElementField(mesh, *current_position,
				    X, *it, _not_ghost);

    Array<Real>::iterator< Matrix<Real> > X_el =
      X.begin(spatial_dimension, nb_nodes_per_element);

    for (UInt el = 0; el < nb_element; ++el, ++X_el, ++eibm) {
      elem.element = el;
      Real el_size    = getFEM().getElementInradius(*X_el, *it);
      Real el_dt      = mat_val[(*eibm)(1)]->getStableTimeStep(el_size, elem);
      min_dt = std::min(min_dt, el_dt);
    }
  }

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real energy = 0.;

  /// call update residual on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    energy += (*mat_it)->getPotentialEnergy();
  }

  /// reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&energy, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return energy;
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
    bool count_node = is_local_node && is_not_pbc_slave_node;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if (count_node)
	mv2 += *vel_val * *vel_val * *mass_val;

      vel_val++;
      mass_val++;
    }
    ekin += mv2;
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&ekin, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return ekin * .5;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy(ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(type);

  Array<Real> vel_on_quad(nb_quadrature_points, spatial_dimension);
  Array<UInt> filter_element(1, 1, index);

  getFEM().interpolateOnQuadraturePoints(*velocity, vel_on_quad,
                                         spatial_dimension,
                                         type, _not_ghost,
                                         &filter_element);

  Array<Real>::iterator< Vector<Real> > vit   = vel_on_quad.begin(spatial_dimension);
  Array<Real>::iterator< Vector<Real> > vend  = vel_on_quad.end(spatial_dimension);

  Vector<Real> rho_v2(nb_quadrature_points);

  Real rho = materials[element_index_by_material(type)(index, 1)]->getRho();

  for (UInt q = 0; vit != vend; ++vit, ++q) {
    rho_v2(q) = rho * vit->dot(*vit);
  }

  AKANTU_DEBUG_OUT();

  return 0.5*getFEM().integrate(rho_v2, type, index);
}


/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getExternalWork() {
  AKANTU_DEBUG_IN();

  Real * velo = velocity->storage();
  Real * forc = force->storage();
  Real * resi = residual->storage();
  bool * boun = boundary->storage();

  Real work = 0.;

  UInt nb_nodes = mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool is_not_pbc_slave_node = !getIsPBCSlaveNode(n);
    bool count_node = is_local_node && is_not_pbc_slave_node;

    for (UInt i = 0; i < spatial_dimension; ++i) {
      if (count_node) {
	if(*boun)
	  work -= *resi * *velo * time_step;
	else
	  work += *forc * *velo * time_step;
      }

      ++velo;
      ++forc;
      ++resi;
      ++boun;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&work, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return work;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id) {
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy();
  } else if (energy_id == "external work"){
    return getExternalWork();
  }

  Real energy = 0.;
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    energy += (*mat_it)->getEnergy(id);
  }

  /// reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&energy, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}
/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id,
				    ElementType & type,
				    UInt index){
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy(type, index);
  }

  std::vector<Material *>::iterator mat_it;
  Vector<UInt> mat = element_index_by_material(type, _not_ghost).begin(2)[index];
  Real energy = materials[mat(1)]->getEnergy(energy_id, type, mat(0));

  AKANTU_DEBUG_OUT();
  return energy;
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesAdded(const Array<UInt> & nodes_list,
				       __attribute__((unused)) const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  if(displacement) displacement->resize(nb_nodes);
  if(mass        ) mass        ->resize(nb_nodes);
  if(velocity    ) velocity    ->resize(nb_nodes);
  if(acceleration) acceleration->resize(nb_nodes);
  if(force       ) force       ->resize(nb_nodes);
  if(residual    ) residual    ->resize(nb_nodes);
  if(boundary    ) boundary    ->resize(nb_nodes);

  if(increment_acceleration) increment_acceleration->resize(nb_nodes);
  if(increment) increment->resize(nb_nodes);

  delete dof_synchronizer;
  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onNodesAdded(nodes_list, event);
  }

  if (method != _explicit_dynamic) {
    delete stiffness_matrix;
    delete jacobian_matrix;
    delete solver;
    SolverOptions  solver_options;
    initImplicit((method == _implicit_dynamic), solver_options);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsAdded(const Array<Element> & element_list,
					  const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  getFEM().initShapeFunctions(_not_ghost);
  getFEM().initShapeFunctions(_ghost);

  Array<Element>::const_iterator<Element> it  = element_list.begin();
  Array<Element>::const_iterator<Element> end = element_list.end();

  /// \todo have rules to choose the correct material
  UInt mat_id = 0;

  UInt * mat_id_vect = NULL;
  try {
    const NewMaterialElementsEvent & event_mat = dynamic_cast<const NewMaterialElementsEvent &>(event);
    mat_id_vect = event_mat.getMaterialList().storage();
  } catch(...) {  }

  for (UInt el = 0; it != end; ++it, ++el) {
    const Element & elem = *it;
    //    element_material(elem.type, elem.ghost_type).push_back(UInt(0));

    if(mat_id_vect) mat_id = mat_id_vect[el];
    Material & mat = *materials[mat_id];

    UInt mat_index = mat.addElement(elem.type, elem.element, elem.ghost_type);
    UInt id[2];
    id[0] = mat_index; id[1] = 0;
    element_index_by_material(elem.type, elem.ghost_type).push_back(id);
  }

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onElementsAdded(element_list, event);
  }

  if(method != _explicit_dynamic) AKANTU_DEBUG_TO_IMPLEMENT();

  assembleMassLumped();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsRemoved(__attribute__((unused)) const Array<Element> & element_list,
					    const ByElementTypeUInt & new_numbering,
					    const RemovedElementsEvent & event) {
  //  MeshUtils::purifyMesh(mesh);

  getFEM().initShapeFunctions(_not_ghost);
  getFEM().initShapeFunctions(_ghost);

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onElementsRemoved(element_list, new_numbering, event);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesRemoved(__attribute__((unused)) const Array<UInt> & element_list,
					 const Array<UInt> & new_numbering,
					 __attribute__((unused)) const RemovedNodesEvent & event) {
  if(displacement) mesh.removeNodesFromArray(*displacement, new_numbering);
  if(mass        ) mesh.removeNodesFromArray(*mass        , new_numbering);
  if(velocity    ) mesh.removeNodesFromArray(*velocity    , new_numbering);
  if(acceleration) mesh.removeNodesFromArray(*acceleration, new_numbering);
  if(force       ) mesh.removeNodesFromArray(*force       , new_numbering);
  if(residual    ) mesh.removeNodesFromArray(*residual    , new_numbering);
  if(boundary    ) mesh.removeNodesFromArray(*boundary    , new_numbering);

  if(increment_acceleration) mesh.removeNodesFromArray(*increment_acceleration, new_numbering);
  if(increment)              mesh.removeNodesFromArray(*increment             , new_numbering);

  delete dof_synchronizer;
  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::addDumpField(const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
#define ADD_FIELD(field, type)						\
  addDumpFieldToDumper(BOOST_PP_STRINGIZE(field),			\
		       new DumperIOHelper::NodalField<type>(*field))

  if(field_id == "displacement")      { ADD_FIELD(displacement, Real); }
  else if(field_id == "mass"        ) { ADD_FIELD(mass        , Real); }
  else if(field_id == "velocity"    ) { ADD_FIELD(velocity    , Real); }
  else if(field_id == "acceleration") { ADD_FIELD(acceleration, Real); }
  else if(field_id == "force"       ) { ADD_FIELD(force       , Real); }
  else if(field_id == "residual"    ) { ADD_FIELD(residual    , Real); }
  else if(field_id == "boundary"    ) { ADD_FIELD(boundary    , bool); }
  else if(field_id == "partitions"  ) {
    addDumpFieldToDumper(field_id,
			 new DumperIOHelper::ElementPartitionField(mesh,
								   spatial_dimension,
								   _not_ghost,
								   _ek_regular));
  }
  else if(field_id == "element_index_by_material") {
    addDumpFieldToDumper(field_id,
			 new DumperIOHelper::ElementalField<UInt>(element_index_by_material,
								  spatial_dimension,
								  _not_ghost,
								  _ek_regular));
  } else {
    addDumpFieldToDumper(field_id,
			 new DumperIOHelper::HomogenizedField<Real,
							      DumperIOHelper::InternalMaterialField>(*this,
												     field_id,
												     spatial_dimension,
												     _not_ghost,
												     _ek_regular));
  }
#undef ADD_FIELD
#endif
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::addDumpFieldArray(const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
#define ADD_FIELD(field, type)						\
  DumperIOHelper::Field * f =						\
    new DumperIOHelper::NodalField<type>(*field);			\
  f->setPadding(3);							\
  addDumpFieldToDumper(BOOST_PP_STRINGIZE(field), f)

  if(field_id == "displacement")      { ADD_FIELD(displacement, Real); }
  else if(field_id == "mass"        ) { ADD_FIELD(mass        , Real); }
  else if(field_id == "velocity"    ) { ADD_FIELD(velocity    , Real); }
  else if(field_id == "acceleration") { ADD_FIELD(acceleration, Real); }
  else if(field_id == "force"       ) { ADD_FIELD(force       , Real); }
  else if(field_id == "residual"    ) { ADD_FIELD(residual    , Real); }
  else {
    typedef DumperIOHelper::HomogenizedField<Real,
					     DumperIOHelper::InternalMaterialField> Field;
    Field * field = new Field(*this,
			      field_id,
			      spatial_dimension,
			      spatial_dimension,
			      _not_ghost,
			      _ek_regular);
    field->setPadding(3);
    addDumpFieldToDumper(field_id, field);
  }
#undef ADD_FIELD
#endif
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::addDumpFieldTensor(const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
  if(field_id == "stress") {
    typedef DumperIOHelper::HomogenizedField<Real,
					     DumperIOHelper::InternalMaterialField,
					     DumperIOHelper::material_stress_field_iterator,
					     DumperIOHelper::AvgHomogenizingFunctor,
					     Matrix> Field;
    Field * field = new Field(*this,
			      field_id,
			      spatial_dimension,
			      spatial_dimension,
			      _not_ghost,
			      _ek_regular);
    field->setPadding(3, 3);
    addDumpFieldToDumper(field_id, field);
  } else if(field_id == "strain") {
    typedef DumperIOHelper::HomogenizedField<Real,
					     DumperIOHelper::InternalMaterialField,
					     DumperIOHelper::material_strain_field_iterator,
					     DumperIOHelper::AvgHomogenizingFunctor,
					     Matrix> Field;
    Field * field = new Field(*this,
			      field_id,
			      spatial_dimension,
			      spatial_dimension,
			      _not_ghost,
			      _ek_regular);
    field->setPadding(3, 3);
    addDumpFieldToDumper(field_id, field);
  } else {
    typedef DumperIOHelper::HomogenizedField<Real,
					     DumperIOHelper::InternalMaterialField,
					     DumperIOHelper::internal_material_field_iterator,
					     DumperIOHelper::AvgHomogenizingFunctor,
					     Matrix> Field;

    Field * field = new Field(*this,
			      field_id,
			      spatial_dimension,
			      spatial_dimension,
			      _not_ghost,
			      _ek_regular);
    field->setPadding(3, 3);
    addDumpFieldToDumper(field_id, field);
  }
#endif
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
  element_index_by_material.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */


__END_AKANTU__
