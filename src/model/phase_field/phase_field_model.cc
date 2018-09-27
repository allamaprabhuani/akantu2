/**
 * @file   phase_field_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Aug 01 2018
 * @date last modification: Wed Aug 01 2018
 *
 * @brief  Implementation of PhaseFieldModel class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "phase_field_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "generalized_trapezoidal.hh"
#include "fe_engine_template.hh"
#include "generalized_trapezoidal.hh"
#include "group_manager_inline_impl.cc"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "parser.hh"
#include "shape_lagrange.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */
namespace akantu {

class Compute
  
  // material functor can be created 
  

  

/* -------------------------------------------------------------------------- */
PhaseFieldModel::PhaseFieldModel(Mesh & mesh, UInt dim, const ID & id,
				 const MemoryID & memory_id)
  : Model(mesh, ModelType::_phase_field_model, dim, id, memory_id),
    damage_gradient("damage_gradient", id),
    damage_on_qpoints("damage_on_qpoints", id),
    strain_history_on_qpoints("strain_history_on_qpoints", id),
    strain_on_qpoints("strain_on_qpoints", id){

  AKANTU_DEBUG_IN();

  this->initDOFManager();
 
  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, _gst_pfm_damage);
    this->registerSynchronizer(synchronizer, _gst_pfm_gradient_damage);
  }

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh, spatial_dimension);

#ifdef AKANTU_USE_IOHELPER
  this->mesh.registerDumper<DumperParaview>("phase_field", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif // AKANTU_USE_IOHELPER

  this->registerParam("l0", l_0, 0., _pat_parsmod, "length scale");
  this->registerParam("gc", g_c, _pat_parsmod, "critical local fracture energy density");
  this->registerParam("lambda", lame_lambda, _pat_parsmod, "lame parameter lambda");
  this->registerParam("mu", lame_mu, _pat_parsmod, "lame paramter mu");

 AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PhaseFieldModel::~PhaseFieldModel() = default;
 

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  damage_on_qpoints.initialize(fem, _nb_component = 1);
  damage_gradient.initialize(fem, _nb_component = spatial_dimension);
  strain_on_qpoints.initialize(fem, _nb_component = spatial_dimension);
  strain_history_on_qpoints.initialize(fem, _nb_component = 1);
  
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  readMaterials();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::readMaterials() {
  auto sect = this->getParserSection();

  if (not std::get<1>(sect)) {
    const auto & section = std::get<0>(sect);
    this->parseSection(section);
  }

  damage_energy_on_qpoints.set(g_c * l_0);
}


  
/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleMatrix(const ID & matrix_id) {

  this->computeStrainHistoryOnQuadPoints(_not_ghost);
  
  if (matrix_id == "K") {
    this->assembleDamageMatrix();
  } else if (matrix_id == "M") { // (and need_to_reassemble_capacity) {
    this->assembleDamageGradMatrix();
  }
}


/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initSolver(TimeStepSolverType time_step_solver_type,
                                 NonLinearSolverType) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->damage, 1, "damage");
  this->allocNodalField(this->external_force, 1, "external_force");
  this->allocNodalField(this->internal_force, 1, "internal_force");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");

  if (!dof_manager.hasDOFs("damage")) {
    dof_manager.registerDOFs("damage", *this->damage, _dst_nodal);
    dof_manager.registerBlockedDOFs("damage", *this->blocked_dofs);
  }

  /*if (time_step_solver_type == _tsst_dynamic ||
      time_step_solver_type == _tsst_dynamic_lumped) {
    this->allocNodalField(this->damage_gradient, spatial_dimension, "damage_gradient");

    if (!dof_manager.hasDOFsDerivatives("damage", 1)) {
      dof_manager.registerDOFsDerivative("damage", 1,
                                         *this->damage_gradient);
    }
    }*/
}



/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleDamageMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new damage matrix");

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }
  this->getDOFManager().clearMatrix("K");

  switch (mesh.getSpatialDimension()) {
  case 1:
    this->assembleDamageMatrix<1>(_not_ghost);
    break;
  case 2:
    this->assembleDamageMatrix<2>(_not_ghost);
    break;
  case 3:
    this->assembleDamageMatrix<3>(_not_ghost);
    break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleDamageGradMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new damage gradient matrix");

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }
  this->getDOFManager().clearMatrix("K");

  switch (mesh.getSpatialDimension()) {
  case 1:
    this->assembleDamageGradMatrix<1>(_not_ghost);
    break;
  case 2:
    this->assembleDamageGradMatrix<2>(_not_ghost);
    break;
  case 3:
    this->assembleDamageGradMatrix<3>(_not_ghost);
    break;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
PhaseFieldModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped", _tsst_dynamic_lumped);
  }
  case _static: {
    return std::make_tuple("static", _tsst_static);
  }
  case _implicit_dynamic: {
    return std::make_tuple("implicit", _tsst_dynamic);
  }
  default:
    return std::make_tuple("unknown", _tsst_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions PhaseFieldModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case _tsst_dynamic_lumped: {
    options.non_linear_solver_type = _nls_lumped;
    options.integration_scheme_type["damage"] = _ist_forward_euler;
    options.solution_type["damage"] = IntegrationScheme::_damage;
    break;
  }
  case _tsst_static: {
    options.non_linear_solver_type = _nls_newton_raphson;
    options.integration_scheme_type["damage"] = _ist_pseudo_time;
    options.solution_type["damage"] = IntegrationScheme::_not_defined;
    break;
  }
  case _tsst_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = _nls_newton_raphson;
      options.integration_scheme_type["damage"] = _ist_forward_euler;
      options.solution_type["damage"] =  IntegrationScheme::_not_defined;
    } else {
      options.non_linear_solver_type = _nls_newton_raphson;
      options.integration_scheme_type["damage"] = _ist_backward_euler;
      options.solution_type["damage"] = IntegrationScheme::_damage;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

  
/* -------------------------------------------------------------------------- */
template<UInt dim>
void PhaseFieldModel::assembleDamageMatrix(const GhostType & ghost_type) {

  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<FEEngineType>();
  ComputeDamageFunctor compute_damage(*this);
  
  for (auto && type : mesh.elementTypes(spatial_dimension, ghost_type)) {

    fem.assembleFieldMatrix(compute_damage, "M", "damage",
			    this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

  
/* -------------------------------------------------------------------------- */
template<UInt dim>
void PhaseFieldModel::assembleDamageGradMatrix(const GhostType & ghost_type) {

  AKANTU_DEBUG_IN();

 
  auto & fem = this->getFEEngine();

  for (auto && type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto bt_d_b = std::make_unique<Array<Real>>(
       nb_element * nb_quadrature_points,
       nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    // damage_energy_on_qpoints = gc*l0 = scalar
    fem.computeBtDB(damage_energy_on_qpoints(type, ghost_type), *bt_d_b, 2, type,
		    ghost_type);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{B}^t * \mathbf{W} *
    /// \mathbf{B}@f$
    auto K_b = std::make_unique<Array<Real>>(
	nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_b");

    fem.integrate(*bt_d_b, *K_b, nb_nodes_per_element * nb_nodes_per_element,
		  type, ghost_type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
       "K", "phasefield", *K_b, type, ghost_type, _symmetric);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PhaseFieldModel::computeStrainHistoryOnQuadPoints(
     const GhostType & ghost_type) {

  Matrix<Real> epsilon_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> epsilon_minus(spatial_dimension, spatial_dimension);
  Matrix<Real> epsilon_dir(spatial_dimension, spatial_dimension);
  Matrix<Real> epsilon_diag_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> epsilon_diag_minus(spatial_dimension, spatial_dimension);

  Vector<Real> eigen_values(spatial_dimension);

  Real trace_plus, trace_minus, phi_plus;
  
  for (auto & type : mesh.elementTypes(spatial_dimension, ghost_type)) {
   
    for (auto  && values :
	   zip(make_view(strain_on_qpoints(type, ghost_type),
			 spatial_dimension, spatial_dimension),
	       make_view(strain_history_on_qpoints(type, ghost_type)))) {

      auto & epsilon     = std::get<0>(values);
      auto & phi_history = std::get<1>(values);
      
      epsilon_plus.clear();
      epsilon_minus.clear();
      epsilon_dir.clear();
      eigen_values.clear();          
      epsilon_diag_plus.clear();
      epsilon_diag_minus.clear();

      epsilon.eig(eigen_values, epsilon_dir);

      for (UInt i=0; i < spatial_dimension; i++) {
	epsilon_diag_plus(i, i)  = std::max(Real(0.), eigen_values(i));
	epsilon_diag_minus(i, i) = std::min(Real(0.), eigen_values(i));
      }

      Matrix<Real> mat_tmp(spatial_dimension, spatial_dimension);
      Matrix<Real> sigma_plus(spatial_dimension, spatial_dimension);
      Matrix<Real> sigma_minus(spatial_dimension, spatial_dimension);

      mat_tmp.mul<false,true>(epsilon_diag_plus, epsilon_dir);
      epsilon_plus.mul<false, false>(epsilon_dir, mat_tmp);
      mat_tmp.mul<false, true>(epsilon_diag_minus, epsilon_dir);
      epsilon_minus.mul<false, true>(epsilon_dir, mat_tmp);

      trace_plus  = std::max(Real(0.), epsilon.trace());
      trace_minus = std::min(Real(0.), epsilon.trace());
      
      for (UInt i=0; i < spatial_dimension; i++) {
	for (UInt j=0; j < spatial_dimension; j++) {
	  sigma_plus(i, j)  = (i==j) * lame_lambda * trace_plus
	                    + 2 * lame_mu * epsilon_plus(i, j);
	  sigma_minus(i, j) = (i==j) * lame_lambda * trace_minus
	                    + 2 * lame_mu * epsilon_minus(i, j);
	}
      }     
      
      phi_plus = 1./2 * sigma_plus.doubleDot(epsilon);

      if (phi_plus > phi_history) {
	phi_history = phi_plus;
      }
     
    }
  }

}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleResidual() {

  AKANTU_DEBUG_IN();

  this->assembleInternalForces();

  this->getDOFManager().assembleToResidual("damage",
					   *this->external_force, 1);
  this->getDOFManager().assembleToResidual("damage",
					   *this->internal_force, 1);


  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  this->internal_force->clear();

  this->synchronize(_gst_pfm_damage);
  auto & fem = this->getFEEngine();

  for (auto ghost_type: ghost_types) {

    computeDrivingForce(ghost_type);
    
    for (auto type: mesh.elementTypes(spatial_dimension, ghost_type)) {
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      auto & driving_force_on_qpoints_vect = driving_force_on_qpoints(type, ghost_type);

      UInt nb_quad_points = driving_force_on_qpoints_vect.size();
      Array<Real> nt_driving_force(nb_quad_points, nb_nodes_per_element);
      fem.computeNtb(driving_force_on_qpoints_vect, nt_driving_force, type, ghost_type);

      UInt nb_elements = mesh.getNbElement(type, ghost_type);
      Array<Real> int_nt_driving_force(nb_elements, nb_nodes_per_element);

      fem.integrate(nt_driving_force, int_nt_driving_force, nb_nodes_per_element, type,
		    ghost_type);

      this->getDOFManager().assembleElementalArrayLocalArray(
	  int_nt_driving_force, *this->internal_force, type, ghost_type, -1);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::computeDrivingForce(const GhostType & ghost_type) {

  for (auto & type: mesh.elementTypes(spatial_dimension, ghost_type)) {

    for (auto && values:
	   zip(make_view(strain_history_on_qpoints(type, ghost_type)),
	       make_view(driving_force_on_qpoints(type, ghost_type)) )) {
      auto & strain_history = std::get<0>(values);
      auto & driving_force  = std::get<1>(values);

      driving_force = 2.0 * strain_history;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt PhaseFieldModel::getNbData(const Array<Element> & elements,
				const SynchronizationTag & tag) const {

  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;

  for (const Element & el : elements) {
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case _gst_material_id: {
    size += elements.size() * sizeof(UInt);
    break;
  }
  case _gst_pfm_damage: {
    size += nb_nodes_per_element * sizeof(Real); // damage
    break;
  }
  case _gst_pfm_gradient_damage: {
    size += getNbIntegrationPoints(elements) * spatial_dimension * sizeof(Real);
    size += nb_nodes_per_element * sizeof(Real);
    break;
  }
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::packData(CommunicationBuffer & buffer,
			       const Array<Element> & elements,
			       const SynchronizationTag & tag) const {

  switch (tag) {
  case _gst_pfm_damage: {
    packNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  case _gst_pfm_gradient_damage: {
    packElementalDataHelper(damage_gradient, buffer, elements, true,
			    getFEEngine());
    packNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }  
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::unpackData(CommunicationBuffer & buffer,
				 const Array<Element> & elements,
				 const SynchronizationTag & tag) {

  switch (tag) {
  case _gst_pfm_damage: {
    unpackNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  case _gst_pfm_gradient_damage: {
    unpackElementalDataHelper(damage_gradient, buffer, elements, true,
			      getFEEngine());
    unpackNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }  
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }
  
}

/* -------------------------------------------------------------------------- */
UInt PhaseFieldModel::getNbData(const Array<UInt> & indexes,
				const SynchronizationTag & tag) const {

  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = indexes.size();

  switch (tag) {
  case _gst_pfm_damage: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::packData(CommunicationBuffer & buffer,
			       const Array<UInt> & indexes,
			       const SynchronizationTag & tag) const {

  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case _gst_pfm_damage: {
      buffer << (*damage)(index);
      break;
    }
    default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::unpackData(CommunicationBuffer & buffer,
				 const Array<UInt> & indexes,
				 const SynchronizationTag & tag) {

  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case _gst_pfm_damage: {
      buffer >> (*damage)(index);
      break;
    }
    default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
    }
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

dumper::Field * PhaseFieldModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag ) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs;

  dumper::Field * field = nullptr;
  field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
  }

/* -------------------------------------------------------------------------- */
dumper::Field * PhaseFieldModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION("Capacity lumped is a nodal field now stored in the DOF manager."
		       "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["damage"] = damage;
  
  dumper::Field * field =
    mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
dumper::Field * PhaseFieldModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const UInt & spatial_dimension,
    const ElementKind & element_kind) {

  dumper::Field * field = nullptr;

  if (field_name == "partitions") 
    field = mesh.createElementalField<UInt, dumper::ElementPartitionField>(
       mesh.getConnectivities(), group_name, this->spatial_dimension,
       element_kind);
  else if (field_name == "damage_gradient") {
    ElementTypeMap<UInt> nb_data_per_elem =
      this->mesh.getNbDataPerElem(damage_gradient, element_kind);

    field = mesh.createElementalField<Real, dumper::InternalMaterialField>(
        damage_gradient, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  } else if (field_name == "conductivity") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(damage_energy_on_qpoints, element_kind);

    field = mesh.createElementalField<Real, dumper::InternalMaterialField>(
	damage_energy_on_qpoints, group_name, this->spatial_dimension,
	element_kind, nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */
dumper::Field * PhaseFieldModel::createElementalField(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const ElementKind & element_kind) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
dumper::Field * PhaseFieldModel::createNodalFieldBool(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
dumper::Field * PhaseFieldModel::createNodalFieldReal(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}
#endif

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::dump(const std::string & dumper_name) {
  mesh.dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::dump(const std::string & dumper_name, UInt step) {
  mesh.dump(dumper_name, step);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::dump(const std::string & dumper_name, Real time,
			   UInt step) {
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::dump() { mesh.dump(); }

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::dump(UInt step) { mesh.dump(step); }

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::dump(Real time, UInt step) { mesh.dump(time, step); }  

/* -------------------------------------------------------------------------- */  
} // akantu
