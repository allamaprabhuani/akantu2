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

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_iohelper_tmpl_homogenizing_field.hh"
#endif

/* -------------------------------------------------------------------------- */
HeatTransferModel::HeatTransferModel(Mesh & mesh,
				     UInt dim,
				     const ID & id,
				     const MemoryID & memory_id) :
  Model(mesh, dim, id, memory_id), Dumpable<DumperParaview>(id),
  integrator(new ForwardEuler()),
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

  this->conductivity_variation = 0.0;
  this->T_ref = 0;

  addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);

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
void HeatTransferModel::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEM().getMesh().getNbNodes();

  std::stringstream sstr_temp;      sstr_temp      << id << ":temperature";
  std::stringstream sstr_temp_rate; sstr_temp_rate << id << ":temperature_rate";
  std::stringstream sstr_inc;       sstr_inc       << id << ":increment";
  std::stringstream sstr_ext_flx;   sstr_ext_flx   << id << ":external_flux";
  std::stringstream sstr_residual;  sstr_residual  << id << ":residual";
  std::stringstream sstr_lump;      sstr_lump      << id << ":lumped";
  std::stringstream sstr_boun;      sstr_boun      << id << ":boundary";

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
  residual->copy(*external_heat_rate);

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
					      1 ,*it,ghost_type);

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
void HeatTransferModel::readMaterials(const std::string & filename) {
  Parser parser;
  parser.open(filename);
  std::string opt_param;
  std::string mat_type = parser.getNextSection("heat", opt_param);

  if (mat_type != "") {
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
    for(UInt i = 0; i < 3; i++)
      for(UInt j = 0; j < 3; j++) {
	if (i< spatial_dimension && j < spatial_dimension){
	  str >> conductivity(i, j);
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
    str >> T_ref;
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
  initArrays();
  temperature->clear();
  temperature_rate->clear();
  external_heat_rate->clear();
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

    Eth += getFEM().integrate(Eth, *it);
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
#define IF_ADD_NODAL_FIELD(field, type, pad)				\
  if(field_id == BOOST_PP_STRINGIZE(field)) {				\
    DumperIOHelper::Field * f =						\
      new DumperIOHelper::NodalField<type>(*field);			\
    if(pad) f->setPadding(3);						\
    addDumpFieldToDumper(BOOST_PP_STRINGIZE(field), f);			\
  } else

#define IF_ADD_ELEM_FIELD(field, type, pad)				\
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
    addDumpFieldToDumper(BOOST_PP_STRINGIZE(field), f);			\
  } else
#endif

void HeatTransferModel::addDumpField(const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
  IF_ADD_NODAL_FIELD(temperature       , Real, false)
  IF_ADD_NODAL_FIELD(temperature_rate  , Real, false)
  IF_ADD_NODAL_FIELD(external_heat_rate, Real, false)
  IF_ADD_NODAL_FIELD(residual          , Real, false)
  IF_ADD_NODAL_FIELD(capacity_lumped   , Real, false)
  IF_ADD_NODAL_FIELD(boundary          , bool, false)
  IF_ADD_ELEM_FIELD (temperature_gradient, Real, false)
  if(field_id == "partitions"  ) {
    addDumpFieldToDumper(field_id,
			 new DumperIOHelper::ElementPartitionField<>(mesh,
                                                                     spatial_dimension,
                                                                     _not_ghost,
                                                                     _ek_regular));
  } else {
    AKANTU_DEBUG_ERROR("Field " << field_id << " does not exists in the model " << id);
  }
#endif
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::addDumpFieldArray(const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
  IF_ADD_ELEM_FIELD (temperature_gradient, Real, false)
  {
    AKANTU_DEBUG_ERROR("Field " << field_id << " does not exists in the model " << id
		       << " or is not dumpable as a vector");
  }
#endif
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::addDumpFieldTensor(const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
  AKANTU_DEBUG_ERROR("Field " << field_id << " does not exists in the model " << id
		       << " or is not dumpable as a Tensor");
#endif
}

__END_AKANTU__
