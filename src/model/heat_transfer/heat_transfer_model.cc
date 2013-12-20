/**
 * @file   heat_transfer_model.cc
 *
 * @author Rui Wang <rui.wang@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Sun May 01 19:14:43 2011
 *
 * @brief  Implementation of HeatTransferModel class
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
#include "aka_common.hh"
#include "heat_transfer_model.hh"
#include "aka_math.hh"
#include "aka_common.hh"
#include "fem_template.hh"
#include "mesh.hh"
#include "static_communicator.hh"
#include "parser.hh"
#include "generalized_trapezoidal.hh"

#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#  include "dumper_iohelper_tmpl_homogenizing_field.hh"
#endif

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

const HeatTransferModelOptions default_heat_transfer_model_options(_explicit_lumped_capacity);

/* -------------------------------------------------------------------------- */
HeatTransferModel::HeatTransferModel(Mesh & mesh,
				     UInt dim,
				     const ID & id,
				     const MemoryID & memory_id) :
  Model(mesh, dim, id, memory_id), Dumpable(), Parsable(_st_heat, id),
  integrator(new ForwardEuler()),
  stiffness_matrix(NULL),
  jacobian_matrix(NULL),
  temperature_gradient    ("temperature_gradient", id),
  temperature_on_qpoints  ("temperature_on_qpoints", id),
  conductivity_on_qpoints ("conductivity_on_qpoints", id),
  k_gradt_on_qpoints      ("k_gradt_on_qpoints", id),
  int_bt_k_gT             ("int_bt_k_gT", id),
  bt_k_gT                 ("bt_k_gT", id),
  conductivity(spatial_dimension, spatial_dimension),
  thermal_energy          ("thermal_energy", id) {
  AKANTU_DEBUG_IN();

  createSynchronizerRegistry(this);

  std::stringstream sstr; sstr << id << ":fem";
  registerFEMObject<MyFEMType>(sstr.str(), mesh,spatial_dimension);

  this->temperature= NULL;
  this->residual = NULL;
  this->boundary = NULL;

#ifdef AKANTU_USE_IOHELPER
  this->registerDumper<DumperParaview>("paraview_all", id, true);
  this->addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif

  this->registerParam("conductivity"          , conductivity              , _pat_parsmod);
  this->registerParam("conductivity_variation", conductivity_variation, 0., _pat_parsmod);
  this->registerParam("temperature_reference" , T_ref                 , 0., _pat_parsmod);
  this->registerParam("capacity"              , capacity                  , _pat_parsmod);
  this->registerParam("density"               , density                   , _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initModel() {
  getFEM().initShapeFunctions(_not_ghost);
  getFEM().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initParallel(MeshPartition * partition,
				     DataAccessor * data_accessor) {
  AKANTU_DEBUG_IN();

  if (data_accessor == NULL) data_accessor = this;
  Synchronizer & synch_parallel = createParallelSynch(partition,data_accessor);

  synch_registry->registerSynchronizer(synch_parallel, _gst_htm_capacity);
  synch_registry->registerSynchronizer(synch_parallel, _gst_htm_temperature);
  synch_registry->registerSynchronizer(synch_parallel, _gst_htm_gradient_temperature);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initPBC() {
  AKANTU_DEBUG_IN();

  Model::initPBC();
  PBCSynchronizer * synch = new PBCSynchronizer(pbc_pair);

  synch_registry->registerSynchronizer(*synch, _gst_htm_capacity);
  synch_registry->registerSynchronizer(*synch, _gst_htm_temperature);
  changeLocalEquationNumberForPBC(pbc_pair,1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEM().getMesh().getNbNodes();

  std::stringstream sstr_temp;      sstr_temp      << Model::id << ":temperature";
  std::stringstream sstr_temp_rate; sstr_temp_rate << Model::id << ":temperature_rate";
  std::stringstream sstr_inc;       sstr_inc       << Model::id << ":increment";
  std::stringstream sstr_ext_flx;   sstr_ext_flx   << Model::id << ":external_flux";
  std::stringstream sstr_residual;  sstr_residual  << Model::id << ":residual";
  std::stringstream sstr_lump;      sstr_lump      << Model::id << ":lumped";
  std::stringstream sstr_boun;      sstr_boun      << Model::id << ":boundary";

  temperature        = &(alloc<Real>(sstr_temp.str(),      nb_nodes, 1, REAL_INIT_VALUE));
  temperature_rate   = &(alloc<Real>(sstr_temp_rate.str(), nb_nodes, 1, REAL_INIT_VALUE));
  increment          = &(alloc<Real>(sstr_inc.str(),       nb_nodes, 1, REAL_INIT_VALUE));
  external_heat_rate = &(alloc<Real>(sstr_ext_flx.str(),   nb_nodes, 1, REAL_INIT_VALUE));
  residual           = &(alloc<Real>(sstr_residual.str(),  nb_nodes, 1, REAL_INIT_VALUE));
  capacity_lumped    = &(alloc<Real>(sstr_lump.str(),      nb_nodes, 1, REAL_INIT_VALUE));
  boundary           = &(alloc<bool>(sstr_boun.str(),      nb_nodes, 1, false));

  Mesh::ConnectivityTypeList::const_iterator it;

  /* -------------------------------------------------------------------------- */
  // byelementtype vectors
  getFEM().getMesh().initByElementTypeArray(temperature_on_qpoints,
					     1,
					     spatial_dimension);

  getFEM().getMesh().initByElementTypeArray(temperature_gradient,
					     spatial_dimension,
					     spatial_dimension);

  getFEM().getMesh().initByElementTypeArray(conductivity_on_qpoints,
					     spatial_dimension*spatial_dimension,
					     spatial_dimension);

  getFEM().getMesh().initByElementTypeArray(k_gradt_on_qpoints,
					     spatial_dimension,
					     spatial_dimension);

  getFEM().getMesh().initByElementTypeArray(bt_k_gT,
					     1,
					     spatial_dimension,true);

  getFEM().getMesh().initByElementTypeArray(int_bt_k_gT,
					     1,
					     spatial_dimension,true);

  getFEM().getMesh().initByElementTypeArray(thermal_energy,
					     1,
					     spatial_dimension);


  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const Mesh::ConnectivityTypeList & type_list =
      getFEM().getMesh().getConnectivityTypeList(gt);

    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
      UInt nb_element = getFEM().getMesh().getNbElement(*it, gt);
      UInt nb_quad_points = this->getFEM().getNbQuadraturePoints(*it, gt) * nb_element;

      temperature_on_qpoints(*it, gt).resize(nb_quad_points);
      temperature_on_qpoints(*it, gt).clear();

      temperature_gradient(*it, gt).resize(nb_quad_points);
      temperature_gradient(*it, gt).clear();

      conductivity_on_qpoints(*it, gt).resize(nb_quad_points);
      conductivity_on_qpoints(*it, gt).clear();

      k_gradt_on_qpoints(*it, gt).resize(nb_quad_points);
      k_gradt_on_qpoints(*it, gt).clear();

      bt_k_gT(*it, gt).resize(nb_quad_points);
      bt_k_gT(*it, gt).clear();

      int_bt_k_gT(*it, gt).resize(nb_element);
      int_bt_k_gT(*it, gt).clear();

      thermal_energy(*it, gt).resize(nb_element);
      thermal_energy(*it, gt).clear();
    }
  }

  /* -------------------------------------------------------------------------- */
  dof_synchronizer = new DOFSynchronizer(getFEM().getMesh(),1);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initSolver(__attribute__((unused)) SolverOptions & options) {
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  delete jacobian_matrix;
  std::stringstream sstr; sstr << Memory::id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(nb_global_nodes, _symmetric,
                                     1, sstr.str(), memory_id);

  jacobian_matrix->buildProfile(mesh, *dof_synchronizer);

  delete stiffness_matrix;
  std::stringstream sstr_sti; sstr_sti << Memory::id << ":stiffness_matrix";
  stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr_sti.str(), memory_id);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << Memory::id << ":solver";
  solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  if(solver)
    solver->initialize(options);
#endif //AKANTU_HAS_SOLVER
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initImplicit(SolverOptions & solver_options) {
  method = _static;
  initSolver(solver_options);
}

/* -------------------------------------------------------------------------- */

HeatTransferModel::~HeatTransferModel()
{
  AKANTU_DEBUG_IN();

  if(dof_synchronizer) delete dof_synchronizer;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  FEM & fem = getFEM();

  const Mesh::ConnectivityTypeList & type_list
    = fem.getMesh().getConnectivityTypeList(ghost_type);

  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it)
  {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_element = getFEM().getMesh().getNbElement(*it,ghost_type);
    UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, ghost_type);

    Array<Real> rho_1 (nb_element * nb_quadrature_points,1, capacity * density);
    fem.assembleFieldLumped(rho_1,1,*capacity_lumped,
			    dof_synchronizer->getLocalDOFEquationNumbers(),
			    *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  capacity_lumped->clear();

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);

  getSynchronizerRegistry().synchronize(_gst_htm_capacity);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::updateResidual() {
  AKANTU_DEBUG_IN();
  /// @f$ r = q_{ext} - q_{int} - C \dot T @f$

  // start synchronization
  synch_registry->asynchronousSynchronize(_gst_htm_temperature);
  // finalize communications
  synch_registry->waitEndSynchronize(_gst_htm_temperature);

  //clear the array
  /// first @f$ r = q_{ext} @f$
  //  residual->clear();
  residual->copy(*external_heat_rate);

  /// then @f$ r -= q_{int} @f$
  // update the not ghost ones

  updateResidual(_not_ghost);
  // update for the received ghosts
  updateResidual(_ghost);

/*  if (method == _explicit_lumped_capacity) {
    this->solveExplicitLumped();
  }*/


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  stiffness_matrix->clear();

  switch(mesh.getSpatialDimension()) {
    case 1: this->assembleStiffnessMatrix<1>(_not_ghost); break;
    case 2: this->assembleStiffnessMatrix<2>(_not_ghost); break;
    case 3: this->assembleStiffnessMatrix<3>(_not_ghost); break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void HeatTransferModel::assembleStiffnessMatrix(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    this->assembleStiffnessMatrix<dim>(*it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */
template <UInt dim>
void HeatTransferModel::assembleStiffnessMatrix(const ElementType & type, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  SparseMatrix & K = *stiffness_matrix;
  const Array<Real> & shapes_derivatives = this->getFEM().getShapesDerivatives(type, ghost_type);

  UInt nb_element                 = mesh.getNbElement(type, ghost_type);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(type, ghost_type);

  Array<Real> & t_gradient = temperature_gradient(type, ghost_type);
  this->getFEM().gradientOnQuadraturePoints(*temperature, t_gradient,
                                             1, type, ghost_type);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = nb_nodes_per_element;

  Array<Real> * bt_d_b = new Array<Real>(nb_element * nb_quadrature_points,
					                               bt_d_b_size * bt_d_b_size,
					                               "B^t*D*B");

  Matrix<Real> Bt_D(nb_nodes_per_element, dim);

  Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_it = shapes_derivatives.begin(dim, nb_nodes_per_element);

  Array<Real>::iterator< Matrix<Real> > Bt_D_B_it = bt_d_b->begin(bt_d_b_size, bt_d_b_size);

  this->computeConductivityOnQuadPoints(ghost_type);
  Array<Real>::iterator< Matrix<Real> > D_it = conductivity_on_qpoints(type, ghost_type).begin(dim, dim);
  Array<Real>::iterator< Matrix<Real> > D_end = conductivity_on_qpoints(type, ghost_type).end(dim, dim);

  for (; D_it != D_end; ++D_it, ++Bt_D_B_it, ++shapes_derivatives_it) {
    Matrix<Real> & D = *D_it;
    const Matrix<Real> & B = *shapes_derivatives_it;
    Matrix<Real> & Bt_D_B = *Bt_D_B_it;

    Bt_D.mul<true, false>(B, D);
    Bt_D_B.mul<false, false>(Bt_D, B);
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> * K_e = new Array<Real>(nb_element,
				      bt_d_b_size * bt_d_b_size,
				      "K_e");

  this->getFEM().integrate(*bt_d_b, *K_e,
                            bt_d_b_size * bt_d_b_size,
                            type, ghost_type);

  delete bt_d_b;

  this->getFEM().assembleMatrix(*K_e, K, 1, type, ghost_type);
  delete K_e;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::solveStatic() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving Ku = f");
  AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
                      "You should first initialize the implicit solver and assemble the stiffness matrix");

  UInt nb_nodes = temperature->getSize();
  UInt nb_degree_of_freedom = temperature->getNbComponent() * nb_nodes;

  jacobian_matrix->copyContent(*stiffness_matrix);
  jacobian_matrix->applyBoundary(*boundary);

  increment->clear();

  solver->setRHS(*residual);
  solver->solve(*increment);

  Real * increment_val = increment->storage();
  Real * temperature_val = temperature->storage();
  bool * boundary_val = boundary->storage();

  for (UInt j = 0; j < nb_degree_of_freedom;
       ++j, ++temperature_val, ++increment_val, ++boundary_val) {
    if (!(*boundary_val)) {
      *temperature_val += *increment_val;
    }
    else {
      *increment_val = 0.0;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeConductivityOnQuadPoints(const GhostType & ghost_type) {
  const Mesh::ConnectivityTypeList & type_list =
    this->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Array<Real> & temperature_interpolated = temperature_on_qpoints(*it, ghost_type);

    //compute the temperature on quadrature points
    this->getFEM().interpolateOnQuadraturePoints(*temperature,
						 temperature_interpolated,
						 1 ,*it,ghost_type);

    Array<Real>::iterator< Matrix<Real> > C_it =
      conductivity_on_qpoints(*it, ghost_type).begin(spatial_dimension, spatial_dimension);
    Array<Real>::iterator< Matrix<Real> > C_end =
      conductivity_on_qpoints(*it, ghost_type).end(spatial_dimension, spatial_dimension);

    Array<Real>::iterator<Real>  T_it = temperature_interpolated.begin();

    for (;C_it != C_end; ++C_it, ++T_it) {
      Matrix<Real> & C = *C_it;
      Real & T = *T_it;
      C = conductivity;

      Matrix<Real> variation(spatial_dimension, spatial_dimension, conductivity_variation * (T - T_ref));
      C += conductivity_variation;
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeKgradT(const GhostType & ghost_type) {

  computeConductivityOnQuadPoints(ghost_type);

  const Mesh::ConnectivityTypeList & type_list =
    this->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    const ElementType & type = *it;
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Array<Real> & gradient = temperature_gradient(*it, ghost_type);
    this->getFEM().gradientOnQuadraturePoints(*temperature,
					      gradient,
					      1 ,*it, ghost_type);

    Array<Real>::iterator< Matrix<Real> > C_it =
      conductivity_on_qpoints(*it, ghost_type).begin(spatial_dimension, spatial_dimension);

    Array<Real>::iterator< Vector<Real> > BT_it = gradient.begin(spatial_dimension);

    Array<Real>::iterator< Vector<Real> > k_BT_it =
      k_gradt_on_qpoints(type, ghost_type).begin(spatial_dimension);
    Array<Real>::iterator< Vector<Real> > k_BT_end =
      k_gradt_on_qpoints(type, ghost_type).end(spatial_dimension);

    for (;k_BT_it != k_BT_end; ++k_BT_it, ++BT_it, ++C_it) {
      Vector<Real> & k_BT = *k_BT_it;
      Vector<Real> & BT = *BT_it;
      Matrix<Real> & C = *C_it;

      k_BT.mul<false>(C, BT);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::updateResidual(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list =
    this->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Array<Real> & shapes_derivatives =
      const_cast<Array<Real> &>(getFEM().getShapesDerivatives(*it,ghost_type));

    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

    // compute k \grad T
    computeKgradT(ghost_type);

    Array<Real>::iterator< Vector<Real> > k_BT_it =
      k_gradt_on_qpoints(*it,ghost_type).begin(spatial_dimension);

    Array<Real>::iterator< Matrix<Real> > B_it =
      shapes_derivatives.begin(spatial_dimension, nb_nodes_per_element);

    Array<Real>::iterator< Vector<Real> > Bt_k_BT_it =
      bt_k_gT(*it,ghost_type).begin(nb_nodes_per_element);

    Array<Real>::iterator< Vector<Real> > Bt_k_BT_end =
      bt_k_gT(*it,ghost_type).end(nb_nodes_per_element);


    for (;Bt_k_BT_it != Bt_k_BT_end; ++Bt_k_BT_it, ++B_it, ++k_BT_it) {
      Vector<Real> & k_BT = *k_BT_it;
      Vector<Real> & Bt_k_BT = *Bt_k_BT_it;
      Matrix<Real> & B = *B_it;

      Bt_k_BT.mul<true>(B, k_BT);
    }

    this->getFEM().integrate(bt_k_gT(*it,ghost_type),
			     int_bt_k_gT(*it,ghost_type),
			     nb_nodes_per_element, *it,ghost_type);

    this->getFEM().assembleArray(int_bt_k_gT(*it,ghost_type), *residual,
				  dof_synchronizer->getLocalDOFEquationNumbers(),
				 1, *it,ghost_type, empty_filter, -1);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::solveExplicitLumped() {
  AKANTU_DEBUG_IN();
  
  /// finally @f$ r -= C \dot T @f$
  // lumped C
  UInt nb_nodes            = temperature_rate->getSize();
  UInt nb_degree_of_freedom = temperature_rate->getNbComponent();

  Real * capacity_val  = capacity_lumped->values;
  Real * temp_rate_val = temperature_rate->values;
  Real * res_val       = residual->values;
  bool * boundary_val  = boundary->values;

  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
    if(!(*boundary_val)) {
      *res_val -= *capacity_val * *temp_rate_val;
    }
    boundary_val++;
    res_val++;
    capacity_val++;
    temp_rate_val++;
  }

#ifndef AKANTU_NDEBUG
  getSynchronizerRegistry().synchronize(akantu::_gst_htm_gradient_temperature);
#endif

  capacity_val      = capacity_lumped->values;
  res_val           = residual->values;
  boundary_val      = boundary->values;
  Real * inc           = increment->values;

  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
    if(!(*boundary_val)) {
      *inc = (*res_val / *capacity_val);
    }
    res_val++;
    boundary_val++;
    inc++;
    capacity_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::explicitPred() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemePred(time_step,
				    *temperature,
				    *temperature_rate,
				    *boundary);

  UInt nb_nodes = temperature->getSize();
  UInt nb_degree_of_freedom = temperature->getNbComponent();

  Real * temp = temperature->values;
  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n, ++temp)
    if(*temp < 0.) *temp = 0.;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::explicitCorr() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemeCorrTempRate(time_step,
					    *temperature,
					    *temperature_rate,
					    *boundary,
					    *increment);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getStableTimeStep()
{
  AKANTU_DEBUG_IN();

  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real conductivitymax = conductivity(0, 0);

  //get the biggest parameter from k11 until k33//
  for(UInt i = 0; i < spatial_dimension; i++)
    for(UInt j = 0; j < spatial_dimension; j++)
      conductivitymax = std::max(conductivity(i, j), conductivitymax);


  const Mesh::ConnectivityTypeList & type_list =
    getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension)
      continue;

    UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(*it);

    Array<Real> coord(0, nb_nodes_per_element*spatial_dimension);
    FEM::extractNodalToElementField(getFEM().getMesh(), getFEM().getMesh().getNodes(),
				    coord, *it, _not_ghost);

    Array<Real>::iterator< Matrix<Real> > el_coord = coord.begin(spatial_dimension, nb_nodes_per_element);
    UInt nb_element           = getFEM().getMesh().getNbElement(*it);

    for (UInt el = 0; el < nb_element; ++el, ++el_coord) {
	el_size    = getFEM().getElementInradius(*el_coord, *it);
	min_el_size = std::min(min_el_size, el_size);
    }

    AKANTU_DEBUG_INFO("The minimum element size : " << min_el_size
		      << " and the max conductivity is : "
		      << conductivitymax);
  }

  Real min_dt = 2 * min_el_size * min_el_size * density
    * capacity / conductivitymax;

  StaticCommunicator::getStaticCommunicator().allReduce(&min_dt, 1, _so_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::readMaterials() {
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
    sub_sect = parser.getSubSections(_st_heat);

  Parser::const_section_iterator it = sub_sect.first;
  const ParserSection & section = *it;

  this->parseSection(section);
}

/* -------------------------------------------------------------------------- */

void HeatTransferModel::initFull(const std::string & material_file, const ModelOptions & options){
  Model::initFull(material_file, options);

  readMaterials();

  const HeatTransferModelOptions & my_options = dynamic_cast<const HeatTransferModelOptions &>(options);

  //initialize the vectors
  initArrays();
  temperature->clear();
  temperature_rate->clear();
  external_heat_rate->clear();

  method = my_options.analysis_method;

  if (method == _static) {
    initImplicit();
  }
}

/* -------------------------------------------------------------------------- */

void HeatTransferModel::initFEMBoundary(bool create_surface) {

  if(create_surface)
    MeshUtils::buildFacets(getFEM().getMesh());

  FEM & fem_boundary = getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::computeThermalEnergyByNode() {
  AKANTU_DEBUG_IN();

  Real ethermal = 0.;

  Array<Real>::iterator< Vector<Real> > heat_rate_it =
    residual->begin(residual->getNbComponent());

  Array<Real>::iterator< Vector<Real> > heat_rate_end =
    residual->end(residual->getNbComponent());


  UInt n = 0;
  for(;heat_rate_it != heat_rate_end; ++heat_rate_it, ++n) {
    Real heat = 0;
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool is_not_pbc_slave_node = !getIsPBCSlaveNode(n);
    bool count_node = is_local_node && is_not_pbc_slave_node;

    Vector<Real> & heat_rate = *heat_rate_it;
    for (UInt i = 0; i < heat_rate.size(); ++i) {
      if (count_node)
	heat += heat_rate[i] * time_step;
    }
    ethermal += heat;
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&ethermal, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return ethermal;
}

/* -------------------------------------------------------------------------- */
template<class iterator>
void HeatTransferModel::getThermalEnergy(iterator Eth,
					 Array<Real>::const_iterator<Real> T_it,
					 Array<Real>::const_iterator<Real> T_end) const {
  for(;T_it != T_end; ++T_it, ++Eth) {
    *Eth = capacity * density * *T_it;
  }
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy(const ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(type);
  Vector<Real> Eth_on_quarature_points(nb_quadrature_points);

  Array<Real>::iterator<Real> T_it  = this->temperature_on_qpoints(type).begin();
  T_it  += index * nb_quadrature_points;

  Array<Real>::iterator<Real> T_end = T_it + nb_quadrature_points;

  getThermalEnergy(Eth_on_quarature_points.storage(), T_it, T_end);

  return getFEM().integrate(Eth_on_quarature_points, type, index);
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy() {
  Real Eth = 0;

  Mesh & mesh = getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  for(; it != last_type; ++it) {
    UInt nb_element = getFEM().getMesh().getNbElement(*it, _not_ghost);
    UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, _not_ghost);
    Array<Real> Eth_per_quad(nb_element * nb_quadrature_points, 1);

    Array<Real>::iterator<Real> T_it  = this->temperature_on_qpoints(*it).begin();
    Array<Real>::iterator<Real> T_end = this->temperature_on_qpoints(*it).end();
    getThermalEnergy(Eth_per_quad.begin(), T_it, T_end);

    Eth += getFEM().integrate(Eth_per_quad, *it);
  }

  return Eth;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id) {
  AKANTU_DEBUG_IN();
  Real energy = 0;

  if("thermal") energy = getThermalEnergy();

  // reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&energy, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & energy_id, const ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();

  Real energy = 0.;

  if("thermal") energy = getThermalEnergy(type, index);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#define IF_ADD_NODAL_FIELD(dumper_name, field, type, pad)		\
  if(field_id == BOOST_PP_STRINGIZE(field)) {				\
    DumperIOHelper::Field * f =						\
      new DumperIOHelper::NodalField<type>(*field);			\
    if(pad) f->setPadding(3);						\
    internalAddDumpFieldToDumper(dumper_name,				\
				 BOOST_PP_STRINGIZE(field), f);		\
  } else

#define IF_ADD_ELEM_FIELD(dumper_name, field, type, pad)		\
  if(field_id == BOOST_PP_STRINGIZE(field)) {				\
    typedef DumperIOHelper::HomogenizedField<type,			\
					     DumperIOHelper::QuadraturePointsField, \
					     DumperIOHelper::quadrature_point_iterator> Field; \
    DumperIOHelper::Field * f =						\
      new Field(getFEM(),						\
		field,							\
		spatial_dimension,					\
		spatial_dimension,					\
		_not_ghost,						\
		_ek_regular);						\
    if(pad) f->setPadding(3);						\
    internalAddDumpFieldToDumper(dumper_name,				\
				 BOOST_PP_STRINGIZE(field), f);		\
  } else
#endif

void HeatTransferModel::addDumpFieldToDumper(const std::string & dumper_name,
					     const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
  IF_ADD_NODAL_FIELD(dumper_name, temperature       , Real, false)
  IF_ADD_NODAL_FIELD(dumper_name, temperature_rate  , Real, false)
  IF_ADD_NODAL_FIELD(dumper_name, external_heat_rate, Real, false)
  IF_ADD_NODAL_FIELD(dumper_name, residual          , Real, false)
  IF_ADD_NODAL_FIELD(dumper_name, capacity_lumped   , Real, false)
  IF_ADD_NODAL_FIELD(dumper_name, boundary          , bool, false)
  IF_ADD_ELEM_FIELD (dumper_name, temperature_gradient, Real, false)
  if(field_id == "partitions"  ) {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
			 new DumperIOHelper::ElementPartitionField<>(mesh,
                                                                     spatial_dimension,
                                                                     _not_ghost,
                                                                     _ek_regular));
  } else {
    AKANTU_DEBUG_ERROR("Field " << field_id << " does not exists in the model " << Model::id);
  }
#endif
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::addDumpFieldVectorToDumper(const std::string & dumper_name, 
						   const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
  IF_ADD_ELEM_FIELD (dumper_name, temperature_gradient, Real, false)
  {
    AKANTU_DEBUG_ERROR("Field " << field_id << " does not exists in the model " << Model::id
		       << " or is not dumpable as a vector");
  }
#endif
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::addDumpFieldTensorToDumper(const std::string & dumper_name,
						   const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
  AKANTU_DEBUG_ERROR("Field " << field_id << " does not exists in the model " << Model::id
		       << " or is not dumpable as a Tensor");
#endif
}

__END_AKANTU__
