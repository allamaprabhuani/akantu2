/**
 * @file   poisson_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Implementation of PoissonModel class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "poisson_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "generalized_trapezoidal.hh"
#include "group_manager_inline_impl.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "parser.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace heat_transfer {
  namespace details {
    class ComputeRhoFunctor {
    public:
      ComputeRhoFunctor(const PoissonModel & model) : model(model){};

      void operator()(Matrix<Real> & rho, const Element & /*unused*/) {
        rho.set(model.getCapacity() * model.getDensity());
      }

    private:
      const PoissonModel & model;
    };
  } // namespace details
} // namespace heat_transfer

/* -------------------------------------------------------------------------- */
PoissonModel::PoissonModel(Mesh & mesh, UInt dim, const ID & id,
			   std::shared_ptr<DOFManager> dof_manager)
  : Model(mesh, ModelType::_poisson_model, std::move(dof_manager), dim, id),
      constitutive_law_index("constitutive law index", id),
      constitutive_law_local_numbering("constitutive law local numbering", id) {
  
  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<FEEngineType>("PoissonModelFEEngine", mesh,
					     Model::spatial_dimension);
   
  this->mesh.registerDumper<DumperParaview>("poisson_model", id, true);
  this->mesh.addDumpMesh(mesh, Model::spatial_dimension, _not_ghost,
			 _ek_regular);

  constitutive_law_selector =
      std::make_shared<DefaultConstitutiveLawSelector>(constitutive_law_index);


  conductivity = Matrix<Real>(this->spatial_dimension, this->spatial_dimension);

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_constitutive_law_id);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_pm_flux);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PoissonModel::~PoissonModel() = default;


/* -------------------------------------------------------------------------- */
void PoissonModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  constitutive_law_index.initialize(mesh, _element_kind = _ek_not_defined,
				    _default_value = UInt(-1), _with_nb_element = true);
  constitutive_law_local_numbering.initialize(mesh, _element_kind = _ek_not_defined,
					      _with_nb_element = true);

  Model::initFullImpl(options);

  // initialize the materials
  if (not this->parser.getLastParsedFile().empty()) {
    this->instantiateConstitutivelaws();
    this->initConstitutiveLaws();
  }

  this->initBC(*this, *dof, *increment, *external_dof_rate);
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType PoissonModel::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions PoissonModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["temperature"] =
        IntegrationSchemeType::_forward_euler;
    options.solution_type["temperature"] = IntegrationScheme::_temperature_rate;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["temperature"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["temperature"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["temperature"] =
          IntegrationSchemeType::_forward_euler;
      options.solution_type["temperature"] =
          IntegrationScheme::_temperature_rate;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["temperature"] =
          IntegrationSchemeType::_backward_euler;
      options.solution_type["temperature"] = IntegrationScheme::_temperature;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
PoissonModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  case _implicit_dynamic: {
    return std::make_tuple("implicit", TimeStepSolverType::_dynamic);
  }
  default:
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }
}


  /* -------------------------------------------------------------------------- */
void PoissonModel::initSolver(TimeStepSolverType time_step_solver_type,
                                   NonLinearSolverType /*unused*/) {
  auto & dof_manager = this->getDOFManager();

  this->allocNodalField(this->dof, 1, "dof");
  this->allocNodalField(this->external_dof_rate, 1, "external_dof_rate");
  this->allocNodalField(this->internal_dof_rate, 1, "internal_dof_rate");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");

  if (!dof_manager.hasDOFs("dof")) {
    dof_manager.registerDOFs("dof", *this->dof, _dst_nodal);
    dof_manager.registerBlockedDOFs("dof", *this->blocked_dofs);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->dof_rate, 1, "dof_rate");

    if (!dof_manager.hasDOFsDerivatives("dof", 1)) {
      dof_manager.registerDOFsDerivative("dof", 1,
                                         *this->dof_rate);
    }
  }
}

  
/* -------------------------------------------------------------------------- */
void PoissonModel::initModel() {
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);

}

  
/* -------------------------------------------------------------------------- */
void PoissonModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  this->assembleInternalDofRate();

  this->getDOFManager().assembleToResidual("dof",
                                           *this->external_dof_rate, 1);
  this->getDOFManager().assembleToResidual("dof",
                                           *this->internal_dof_rate, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MatrixType PoissonModel::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K" or matrix_id == "M") {
    return _symmetric;
  }

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void PoissonModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacity();
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacityLumped();
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::beforeSolveStep() {
  for (auto & law : constitutive_laws) {
    law->beforeSolveStep();
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::afterSolveStep(bool converged) {
  for (auto & law : constitutive_laws) {
    law->afterSolveStep(converged);
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::predictor() { ++temperature_release; }
 
/* -------------------------------------------------------------------------- */
void PoissonModel::corrector() { ++displacement_release; }

/* -------------------------------------------------------------------------- */
void PoissonModel::assembleInternalDofRate() {
  AKANTU_DEBUG_IN();

  this->internal_dof_rate->zero();

  // compute the fluxes of local elements
  AKANTU_DEBUG_INFO("Compute local fluxes");
  for (auto & law : constitutive_laws) {
    law->computeAllFluxes(_not_ghost);
  }

  
  // communicate the stresses
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  this->asynchronousSynchronize(SynchronizationTag::_pm_flux);
  
  // assemble the rate due to local fluxes
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for (auto & law : constitutive_laws) {
    law->assembleInternalDofRate(_not_ghost);
  }

  // finalize communications
  AKANTU_DEBUG_INFO("Wait distant stresses");
  this->waitEndSynchronize(SynchronizationTag::_pm_flux);

    // assemble the fluxes due to ghost elements
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for (auto & law : constitutive_laws) {
    law->assembleInternalDofRate(_ghost);
  }

 
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PoissonModel::assembleStiffnessMatrix(bool need_to_reassemble) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  // Check if constitutive laws need to recompute the matrix
  for (auto & law : constitutive_laws) {
    need_to_reassemble |= law->hasMatrixChanged("K");
  }

  if (need_to_reassemble) {
    this->getDOFManager().getMatrix("K").zero();

    // call compute stiffness matrix on each local elements
    for (auto & law : constitutive_laws) {
      law->assembleStiffnessMatrix(_not_ghost);
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PoissonModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  if (not need_to_reassemble_capacity_lumped) {
    return;
  }

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().zeroLumpedMatrix("M");

  // call compute capacity lumped equivalent matrix on each local elements
  for (auto & law : constitutive_laws) {
      law->assembleCapacityLumped(_not_ghost);
      law->assembleCapacityLumped(_ghost);
  }
  
  need_to_reassemble_capacity_lumped = false;

  AKANTU_DEBUG_OUT();
}
  
  
/* -------------------------------------------------------------------------- */
FEEngine & PoissonModel::getFEEngineBoundary(const ID & name) {
  return aka::as_type<FEEngine>(getFEEngineClassBoundary<FEEngineType>(name));
}


/* -------------------------------------------------------------------------- */
void PoissonModel::assembleCapacityLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<FEEngineType>();
  heat_transfer::details::ComputeRhoFunctor compute_rho(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldLumped(compute_rho, "M", "temperature",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void PoissonModel::assembleConductivityMatrix() {
  AKANTU_DEBUG_IN();

  this->computeConductivityOnQuadPoints(_not_ghost);

  if (conductivity_release[_not_ghost] == conductivity_matrix_release) {
    return;
  }

  AKANTU_DEBUG_ASSERT(this->getDOFManager().hasMatrix("K"),
                      "The K matrix has not been initialized yet.");

  this->getDOFManager().zeroMatrix("K");

  auto & fem = this->getFEEngine();

  for (auto && type : mesh.elementTypes(spatial_dimension)) {
    auto nb_element = mesh.getNbElement(type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type);

    auto bt_d_b = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    fem.computeBtDB(conductivity_on_qpoints(type), *bt_d_b, 2, type);

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    auto K_e = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_e");

    fem.integrate(*bt_d_b, *K_e, nb_nodes_per_element * nb_nodes_per_element,
                  type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "temperature", *K_e, type, _not_ghost, _symmetric);
  }

  conductivity_matrix_release = conductivity_release[_not_ghost];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PoissonModel::computeConductivityOnQuadPoints(GhostType ghost_type) {
  // if already computed once check if need to compute
  if (not initial_conductivity[ghost_type]) {
    // if temperature did not change, conductivity will not vary
    if (temperature_release == conductivity_release[ghost_type]) {
      return;
    }
    // if conductivity_variation is 0 no need to recompute
    if (conductivity_variation == 0.) {
      return;
    }
  }

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & temperature_interpolated = temperature_on_qpoints(type, ghost_type);

    // compute the temperature on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        *temperature, temperature_interpolated, 1, type, ghost_type);

    auto & cond = conductivity_on_qpoints(type, ghost_type);
    for (auto && tuple :
         zip(make_view(cond, spatial_dimension, spatial_dimension),
             temperature_interpolated)) {
      auto & C = std::get<0>(tuple);
      auto & T = std::get<1>(tuple);
      C = conductivity;

      Matrix<Real> variation(spatial_dimension, spatial_dimension,
                             conductivity_variation * (T - T_ref));
      // @TODO: Guillaume are you sure ? why due you compute variation then ?
      C += conductivity_variation;
    }
  }

  conductivity_release[ghost_type] = temperature_release;
  initial_conductivity[ghost_type] = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PoissonModel::computeKgradT(GhostType ghost_type) {
  computeConductivityOnQuadPoints(ghost_type);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & gradient = temperature_gradient(type, ghost_type);
    this->getFEEngine().gradientOnIntegrationPoints(*temperature, gradient, 1,
                                                    type, ghost_type);

    for (auto && values :
         zip(make_view(conductivity_on_qpoints(type, ghost_type),
                       spatial_dimension, spatial_dimension),
             make_view(gradient, spatial_dimension),
             make_view(k_gradt_on_qpoints(type, ghost_type),
                       spatial_dimension))) {
      const auto & C = std::get<0>(values);
      const auto & BT = std::get<1>(values);
      auto & k_BT = std::get<2>(values);

      k_BT.mul<false>(C, BT);
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
Real PoissonModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real conductivitymax = conductivity(0, 0);

  // get the biggest parameter from k11 until k33//
  for (UInt i = 0; i < spatial_dimension; i++) {
    for (UInt j = 0; j < spatial_dimension; j++) {
      conductivitymax = std::max(conductivity(i, j), conductivitymax);
    }
  }
  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {

    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

    Array<Real> coord(0, nb_nodes_per_element * spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), coord, type,
                                         _not_ghost);

    auto el_coord = coord.begin(spatial_dimension, nb_nodes_per_element);
    UInt nb_element = mesh.getNbElement(type);

    for (UInt el = 0; el < nb_element; ++el, ++el_coord) {
      el_size = getFEEngine().getElementInradius(*el_coord, type);
      min_el_size = std::min(min_el_size, el_size);
    }

    AKANTU_DEBUG_INFO("The minimum element size : "
                      << min_el_size
                      << " and the max conductivity is : " << conductivitymax);
  }

  Real min_dt = 2. * min_el_size * min_el_size / 4. * density * capacity /
                conductivitymax;

  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}
/* -------------------------------------------------------------------------- */

void PoissonModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

  this->mesh.getDumper("heat_transfer").setTimeStep(time_step);
}


/* -------------------------------------------------------------------------- */
void PoissonModel::assembleCapacity() {
  AKANTU_DEBUG_IN();
  auto ghost_type = _not_ghost;

  this->getDOFManager().zeroMatrix("M");

  auto & fem = getFEEngineClass<FEEngineType>();

  heat_transfer::details::ComputeRhoFunctor rho_functor(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldMatrix(rho_functor, "M", "temperature",
                            this->getDOFManager(), type, ghost_type);
  }

  need_to_reassemble_capacity = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PoissonModel::computeRho(Array<Real> & rho, ElementType type,
			      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = this->getFEEngine();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  rho.resize(nb_element * nb_quadrature_points);
  rho.set(this->capacity);

  // Real * rho_1_val = rho.storage();
  // /// compute @f$ rho @f$ for each nodes of each element
  // for (UInt el = 0; el < nb_element; ++el) {
  //   for (UInt n = 0; n < nb_quadrature_points; ++n) {
  //     *rho_1_val++ = this->capacity;
  //   }
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real PoissonModel::computeThermalEnergyByNode() {
  AKANTU_DEBUG_IN();

  Real ethermal = 0.;

  for (auto && pair : enumerate(make_view(
           *internal_heat_rate, internal_heat_rate->getNbComponent()))) {
    auto n = std::get<0>(pair);
    auto & heat_rate = std::get<1>(pair);

    Real heat = 0.;
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool count_node = is_local_node;

    for (UInt i = 0; i < heat_rate.size(); ++i) {
      if (count_node) {
        heat += heat_rate[i] * time_step;
      }
    }
    ethermal += heat;
  }

  mesh.getCommunicator().allReduce(ethermal, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return ethermal;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PoissonModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  auto field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PoissonModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION(
        "Capacity lumped is a nodal field now stored in the DOF manager."
        "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["dof"] = dof.get();
  real_nodal_fields["dof_rate"] = dof_rate.get();
  real_nodal_fields["external_dof_rate"] = external_dof_rate.get();
  real_nodal_fields["internal_dof_rate"] = internal_dof_rate.get();
  real_nodal_fields["increment"] = increment.get();

  std::shared_ptr<dumpers::Field> field =
      mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PoissonModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*padding_flag*/, UInt /*spatial_dimension*/,
    ElementKind element_kind) {

  std::shared_ptr<dumpers::Field> field;

  if (field_name == "partitions") {
    field = mesh.createElementalField<UInt, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  } else if (field_name == "temperature_gradient") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(temperature_gradient);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        temperature_gradient, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  } else if (field_name == "conductivity") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(conductivity_on_qpoints);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        conductivity_on_qpoints, group_name, this->spatial_dimension,
        element_kind, nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
inline UInt PoissonModel::getNbData(const Array<UInt> & indexes,
                                         const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = indexes.size();

  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::packData(CommunicationBuffer & buffer,
                                        const Array<UInt> & indexes,
                                        const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_htm_temperature: {
      buffer << (*temperature)(index);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::unpackData(CommunicationBuffer & buffer,
                                          const Array<UInt> & indexes,
                                          const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_htm_temperature: {
      buffer >> (*temperature)(index);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt PoissonModel::getNbData(const Array<Element> & elements,
                                         const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;
  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  case SynchronizationTag::_htm_gradient_temperature: {
    // temperature gradient
    size += getNbIntegrationPoints(elements) * spatial_dimension * sizeof(Real);
    size += nb_nodes_per_element * sizeof(Real); // nodal temperatures
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::packData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_htm_gradient_temperature: {
    packElementalDataHelper(temperature_gradient, buffer, elements, true,
                            getFEEngine());
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::unpackData(CommunicationBuffer & buffer,
                                          const Array<Element> & elements,
                                          const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_htm_gradient_temperature: {
    unpackElementalDataHelper(temperature_gradient, buffer, elements, true,
                              getFEEngine());
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);

    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
