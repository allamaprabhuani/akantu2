/**
 * @file   heat_transfer_model.cc
 *
 * @author Rui Wang <rui.wang@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
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

// #include "sparse_matrix.hh"
// #include "solver.hh"
// #ifdef AKANTU_USE_MUMPS
// #include "solver_mumps.hh"
// #endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
HeatTransferModel::HeatTransferModel(Mesh & mesh,
                      UInt dim,
                      const ID & id,
                      const MemoryID & memory_id) :
  Model(id, memory_id),
  integrator(new ForwardEuler()),
  spatial_dimension(dim),
  temperature_gradient    ("temperature_gradient", id),
  temperature_on_qpoints  ("temperature_on_qpoints", id),
  conductivity_on_qpoints ("conductivity_on_qpoints", id),
  k_gradt_on_qpoints      ("k_gradt_on_qpoints", id),
  int_bt_k_gT             ("int_bt_k_gT", id),
  bt_k_gT                 ("bt_k_gT", id)
{
  AKANTU_DEBUG_IN();

  createSynchronizerRegistry(this);

  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();

  std::stringstream sstr; sstr << id << ":fem";
  registerFEMObject<MyFEMType>(sstr.str(), mesh,spatial_dimension);

  this->temperature= NULL;

  this->residual = NULL;
  this->boundary = NULL;
  this->conductivity_variation = 0.0;
  this->t_ref = 0;

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
  changeLocalEquationNumberforPBC(pbc_pair,1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initVectors() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEM().getMesh().getNbNodes();

  std::stringstream sstr_temp;      sstr_temp     << id << ":temperature";
  std::stringstream sstr_temp_rate; sstr_temp     << id << ":temperature_rate";
  std::stringstream sstr_inc;       sstr_temp_rate<< id << ":increment";
  std::stringstream sstr_ext_flx;   sstr_ext_flx  << id << ":external_flux";
  std::stringstream sstr_residual;  sstr_residual << id << ":residual";
  std::stringstream sstr_lump;      sstr_lump     << id << ":lumped";
  std::stringstream sstr_boun;      sstr_boun     << id << ":boundary";

  temperature      = &(alloc<Real>(sstr_temp.str(),      nb_nodes, 1, REAL_INIT_VALUE));
  temperature_rate = &(alloc<Real>(sstr_temp_rate.str(), nb_nodes, 1, REAL_INIT_VALUE));
  increment        = &(alloc<Real>(sstr_inc.str(),       nb_nodes, 1, REAL_INIT_VALUE));
  external_flux    = &(alloc<Real>(sstr_ext_flx.str(),   nb_nodes, 1, REAL_INIT_VALUE));
  residual         = &(alloc<Real>(sstr_residual.str(),  nb_nodes, 1, REAL_INIT_VALUE));
  capacity_lumped  = &(alloc<Real>(sstr_lump.str(),      nb_nodes, 1, REAL_INIT_VALUE));
  boundary         = &(alloc<bool>(sstr_boun.str(),      nb_nodes, 1, false));

  Mesh::ConnectivityTypeList::const_iterator it;

  /* -------------------------------------------------------------------------- */
  // byelementtype vectors

  getFEM().getMesh().initByElementTypeVector(temperature_on_qpoints,
                                             1,
                                             spatial_dimension);

  getFEM().getMesh().initByElementTypeVector(temperature_gradient,
                                             spatial_dimension,
                                             spatial_dimension);

  getFEM().getMesh().initByElementTypeVector(conductivity_on_qpoints,
                                             spatial_dimension*spatial_dimension,
                                             spatial_dimension);

  getFEM().getMesh().initByElementTypeVector(k_gradt_on_qpoints,
                                             spatial_dimension,
                                             spatial_dimension);

  getFEM().getMesh().initByElementTypeVector(bt_k_gT,
                                             1,
                                             spatial_dimension,true);

  getFEM().getMesh().initByElementTypeVector(int_bt_k_gT,
                                             1,
                                             spatial_dimension,true);


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
    }
  }

  /* -------------------------------------------------------------------------- */
  dof_synchronizer = new DOFSynchronizer(getFEM().getMesh(),1);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

HeatTransferModel::~HeatTransferModel()
{
  AKANTU_DEBUG_IN();

  delete [] conductivity;
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

    Vector<Real> rho_1 (nb_element * nb_quadrature_points,1, capacity * density);
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

  getSynchronizerRegistry().synchronize(akantu::_gst_htm_capacity);

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
  residual->copy(*external_flux);

  /// then @f$ r -= q_{int} @f$
  // update the not ghost ones

  updateResidual(_not_ghost);
  // update for the received ghosts
  updateResidual(_ghost);

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

  solveExplicitLumped();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeConductivityOnQuadPoints(const GhostType & ghost_type) {
  const Mesh::ConnectivityTypeList & type_list =
    this->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    //compute the temperature on quadrature points
    this->getFEM().interpolateOnQuadraturePoints(*temperature,
                                                 temperature_on_qpoints(*it, ghost_type),
                                                 1 ,*it,ghost_type);

    Vector<Real>::iterator<types::RMatrix> c_iterator =
      conductivity_on_qpoints(*it,ghost_type).begin(spatial_dimension,spatial_dimension);

    Real * T_val = temperature_on_qpoints(*it, ghost_type).values;

    UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, ghost_type);
    UInt nb_element = getFEM().getMesh().getNbElement(*it,ghost_type);
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quadrature_points; ++q) {
        Real * c_onquad_val = c_iterator->storage();
        for (UInt i = 0; i < spatial_dimension*spatial_dimension; ++i) {
          //      if (i == 0 && *T_val - t_ref > 0)
          // AKANTU_DEBUG_INFO("conductivity is " << conductivity[i]
          //                << " variation is " << conductivity_variation
          //                << " temp increase is " << *T_val - t_ref
          //                << " new conductivity is " << conductivity[i] + conductivity_variation* (*T_val - t_ref));
          c_onquad_val[i] = conductivity[i] + conductivity_variation* (*T_val - t_ref);
        }
        ++T_val;
        ++c_iterator;
      }
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
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_element = getFEM().getMesh().getNbElement(*it,ghost_type);
    UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, ghost_type);


    Vector<Real> & gradient = temperature_gradient(*it, ghost_type);
    gradient.resize(nb_element * nb_quadrature_points);

    this->getFEM().gradientOnQuadraturePoints(*temperature,
                                              gradient,
                                              1 ,*it,ghost_type);

    Vector<Real>::iterator<types::RMatrix> c_iterator =
      conductivity_on_qpoints(*it,ghost_type).begin(spatial_dimension,spatial_dimension);

    Vector<Real>::iterator<types::RVector> gT_iterator =
      temperature_gradient(*it,ghost_type).begin(spatial_dimension);

    Vector<Real>::iterator<types::RVector> k_gT_iterator =
      k_gradt_on_qpoints(*it,ghost_type).begin(spatial_dimension);

    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt i = 0; i < nb_quadrature_points; ++i){
        Math::matrixt_vector(spatial_dimension, spatial_dimension,
                             c_iterator->storage(),
                             gT_iterator->storage(),
                             k_gT_iterator->storage());
        ++c_iterator;
        ++gT_iterator;
        ++k_gT_iterator;
      }
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

    Vector<Real> & shapes_derivatives =
      const_cast<Vector<Real> &>(getFEM().getShapesDerivatives(*it,ghost_type));

    UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, ghost_type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt nb_element = getFEM().getMesh().getNbElement(*it,ghost_type);

    // compute k \grad T
    computeKgradT(ghost_type);

    Vector<Real>::iterator<types::RVector> k_gT_iterator =
      k_gradt_on_qpoints(*it,ghost_type).begin(spatial_dimension);

    Vector<Real>::iterator<types::RMatrix> shapesd_iterator =
      shapes_derivatives.begin(nb_nodes_per_element,spatial_dimension);

    Vector<Real>::iterator<types::RVector> bt_k_gT_iterator =
      bt_k_gT(*it,ghost_type).begin(nb_nodes_per_element);


    // UInt bt_k_gT_size = nb_nodes_per_element;
    // Vector<Real> * bt_k_gT =
    //   new Vector <Real> (nb_quadrature_points*nb_element, bt_k_gT_size);
    // Real * bt_k_gT_val = bt_k_gT->values;

    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt i = 0; i < nb_quadrature_points; ++i) {
        Math::matrix_vector(nb_nodes_per_element, spatial_dimension,
                            shapesd_iterator->storage(),
                            k_gT_iterator->storage(),
                            bt_k_gT_iterator->storage());


        ++shapesd_iterator;
        ++k_gT_iterator;
        ++bt_k_gT_iterator;
      }
    }

    this->getFEM().integrate(bt_k_gT(*it,ghost_type),
                             int_bt_k_gT(*it,ghost_type),
                             nb_nodes_per_element, *it,ghost_type);

    this->getFEM().assembleVector(int_bt_k_gT(*it,ghost_type), *residual,
                                  dof_synchronizer->getLocalDOFEquationNumbers(),
                                  1, *it,ghost_type,NULL,-1);
    //delete q_e;

  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::solveExplicitLumped() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = increment->getSize();
  UInt nb_degree_of_freedom = increment->getNbComponent();

  Real * capa_val      = capacity_lumped->values;
  Real * res_val       = residual->values;
  bool * boundary_val  = boundary->values;
  Real * inc           = increment->values;

  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
    if(!(*boundary_val)) {
      *inc = (*res_val / *capa_val);
    }
    res_val++;
    boundary_val++;
    inc++;
    capa_val++;
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

  Real conductivitymax = -std::numeric_limits<Real>::max();
  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real * coord = getFEM().getMesh().getNodes().values;

  const Mesh::ConnectivityTypeList & type_list =
    getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension)
      continue;

    UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(*it);
    UInt nb_element           = getFEM().getMesh().getNbElement(*it);

    UInt * conn = getFEM().getMesh().getConnectivity(*it, _not_ghost).values;
    Real * u = new Real[nb_nodes_per_element*spatial_dimension];

    for (UInt el = 0; el < nb_element; ++el) {
        UInt el_offset  = el * nb_nodes_per_element;
        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
          UInt offset_conn = conn[el_offset + n] * spatial_dimension;
          memcpy(u + n * spatial_dimension,
                 coord + offset_conn,
                 spatial_dimension * sizeof(Real));
        }

        el_size    = getFEM().getElementInradius(u, *it);
        min_el_size = std::min(min_el_size, el_size);
    }

    AKANTU_DEBUG_INFO("The minimum element size : " << min_el_size
                      << " and the max conductivity is : "
                      << conductivitymax);
    delete [] u;
  }

  //get the biggest parameter from k11 until k33//
  for(UInt i = 0; i < spatial_dimension * spatial_dimension; i++) {
    if(conductivity[i] > conductivitymax)
      conductivitymax=conductivity[i];
  }

  Real min_dt = 2 * min_el_size * min_el_size * density
    * capacity/conductivitymax;

  StaticCommunicator::getStaticCommunicator().allReduce(&min_dt, 1, _so_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::readMaterials(const std::string & filename) {
  Parser parser;
  parser.open(filename);
  std::string opt_param;
  std::string mat_type = parser.getNextSection("heat", opt_param);

  if (mat_type != ""){
    parser.readSection<HeatTransferModel>(*this);
  }
  else
    AKANTU_DEBUG_ERROR("did not find any section with material info "
                       <<"for heat conduction");
}

/* -------------------------------------------------------------------------- */
bool HeatTransferModel::setParam(const std::string & key,
                                 const std::string & value) {
  std::stringstream str(value);
  if (key == "conductivity") {
    conductivity = new Real[spatial_dimension * spatial_dimension];
    for(UInt i = 0; i < 3; i++)
      for(UInt j = 0; j < 3; j++) {
        if (i< spatial_dimension && j < spatial_dimension){
          str >> conductivity[i*spatial_dimension+j];
          AKANTU_DEBUG_INFO("Conductivity(" << i << "," << j << ") = "
                            << conductivity[i*spatial_dimension+j]);
        } else {
          Real tmp;
          str >> tmp;
        }
      }
  }
  else if (key == "conductivity_variation") {
    str >> conductivity_variation;
  }
  else if (key == "temperature_reference") {
    str >> t_ref;
  }
  else if (key == "capacity"){
    str >> capacity;
    AKANTU_DEBUG_INFO("The capacity of the material is:" << capacity);
  }
  else if (key == "density"){
    str >> density;
    AKANTU_DEBUG_INFO("The density of the material is:" << density);
  } else {
    return false;
  }

  return true;
}

/* -------------------------------------------------------------------------- */

void HeatTransferModel::initFull(const std::string & material_file){
  readMaterials(material_file);
  //model initialization
  initModel();
  //initialize the vectors
  initVectors();
  temperature->clear();
  temperature_rate->clear();
  external_flux->clear();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initFEMBoundary(bool create_surface) {

  if(create_surface)
    MeshUtils::buildFacets(getFEM().getMesh());

  FEM & fem_boundary = getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();
}


__END_AKANTU__
