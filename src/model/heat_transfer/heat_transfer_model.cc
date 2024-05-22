/**
 * @file   heat_transfer_model.cc
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
 * @brief  Implementation of HeatTransferModel class
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
#include "heat_transfer_model.hh"
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
      ComputeRhoFunctor(const HeatTransferModel & model) : model(model){};

      void operator()(Matrix<Real> & rho, const Element & el) {
        auto capacity_it =
            model.getCapacityArray(el.type, el.ghost_type).begin();
        auto density_it = model.getDensityArray(el.type, el.ghost_type).begin();

        rho.set(capacity_it[el.element] * density_it[el.element]);
      }

    private:
      const HeatTransferModel & model;
    };
  } // namespace details
} // namespace heat_transfer

/* -------------------------------------------------------------------------- */
HeatTransferModel::HeatTransferModel(Mesh & mesh, UInt dim, const ID & id,
                                     ModelType model_type,
                                     std::shared_ptr<DOFManager> dof_manager)
    : Model(mesh, model_type, dof_manager, dim, id),
      temperature_gradient("temperature_gradient", id),
      temperature_on_qpoints("temperature_on_qpoints", id),
      temperature_rate_on_qpoints("temperature_rate_on_qpoints", id),
      conductivity_on_qpoints("conductivity_on_qpoints", id),
      k_gradt_on_qpoints("k_gradt_on_qpoints", id),
      density_array("density_array", id), capacity_array("capacity_array", id),
      T_ref_array("temperature_reference_array", id),
      initial_conductivity_array("initial_conductivity_array", id),
      conductivity_variation_array("conductivity_variation_array", id) {
  AKANTU_DEBUG_IN();

  conductivity = Matrix<Real>(this->spatial_dimension, this->spatial_dimension);

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_temperature);
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_temperature_on_qpoints);
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_gradient_temperature);
  }

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh, spatial_dimension);

  this->mesh.registerDumper<DumperParaview>("heat_transfer", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);

  this->registerParam("conductivity", conductivity, _pat_parsmod);
  this->registerParam("conductivity_variation", conductivity_variation, 0.,
                      _pat_parsmod);
  this->registerParam("temperature_reference", T_ref, 0., _pat_parsmod);
  this->registerParam("capacity", capacity, _pat_parsmod);
  this->registerParam("density", density, _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  temperature_on_qpoints.initialize(fem, _nb_component = 1);
  temperature_rate_on_qpoints.initialize(fem, _nb_component = 1);
  temperature_gradient.initialize(fem, _nb_component = spatial_dimension);
  conductivity_on_qpoints.initialize(fem, _nb_component = spatial_dimension *
                                                          spatial_dimension);
  k_gradt_on_qpoints.initialize(fem, _nb_component = spatial_dimension);
  capacity_array.initialize(fem, _nb_component = 1);
  density_array.initialize(fem, _nb_component = 1);
  initial_conductivity_array.initialize(fem, _nb_component = spatial_dimension *
                                                             spatial_dimension);
  conductivity_variation_array.initialize(fem, _nb_component = 1);
  T_ref_array.initialize(fem, _nb_component = 1);

  this->initBC(*this, *temperature, *increment, *external_heat_rate);
}

/* -------------------------------------------------------------------------- */
FEEngine & HeatTransferModel::getFEEngineBoundary(const ID & name) {
  return aka::as_type<FEEngine>(getFEEngineClassBoundary<FEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
HeatTransferModel::~HeatTransferModel() = default;

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped(GhostType ghost_type) {
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
MatrixType HeatTransferModel::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K" or matrix_id == "M") {
    return _symmetric;
  }

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleConductivityMatrix();
  } else if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacity();
  }
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacityLumped();
  }
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  this->assembleInternalHeatRate();

  this->getDOFManager().assembleToResidual("temperature",
                                           *this->external_heat_rate, 1);
  this->getDOFManager().assembleToResidual("temperature",
                                           *this->internal_heat_rate, -1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::predictor() { ++temperature_release; }

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().zeroLumpedMatrix("M");

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);

  need_to_reassemble_capacity_lumped = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initSolver(TimeStepSolverType time_step_solver_type,
                                   NonLinearSolverType /*unused*/) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->temperature, 1, "temperature");
  this->allocNodalField(this->external_heat_rate, 1, "external_heat_rate");
  this->allocNodalField(this->internal_heat_rate, 1, "internal_heat_rate");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");

  if (!dof_manager.hasDOFs("temperature")) {
    dof_manager.registerDOFs("temperature", *this->temperature, _dst_nodal);
    dof_manager.registerBlockedDOFs("temperature", *this->blocked_dofs);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->temperature_rate, 1, "temperature_rate");

    if (!dof_manager.hasDOFsDerivatives("temperature", 1)) {
      dof_manager.registerDOFsDerivative("temperature", 1,
                                         *this->temperature_rate);
    }
  }
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
HeatTransferModel::getDefaultSolverID(const AnalysisMethod & method) {
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
ModelSolverOptions HeatTransferModel::getDefaultSolverOptions(
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
void HeatTransferModel::assembleConductivityMatrix() {
  AKANTU_DEBUG_IN();

  this->computeConductivityOnQuadPoints(_not_ghost);

  if (conductivity_release[_not_ghost] == conductivity_matrix_release) {
    return;
  }

  AKANTU_DEBUG_ASSERT(this->getDOFManager().hasMatrix("K"),
                      "The K matrix has not been initialized yet.");

  this->getDOFManager().zeroMatrix("K");

  auto & fem = this->getFEEngine();

  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
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
        "K", "temperature", *K_e, type, _symmetric);
  }

  conductivity_matrix_release = conductivity_release[_not_ghost];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeConductivityOnQuadPoints(GhostType ghost_type) {
  // if already computed once check if need to compute
  if (not initial_conductivity[ghost_type]) {
    // if temperature did not change, conductivity will not vary
    if (temperature_release == conductivity_release[ghost_type]) {
      return;
    }
    // // if conductivity_variation is 0 no need to recompute
    if (conductivity_variation == 0.) {
      return;
    }
  }

  computeTempOnQpoints(_not_ghost);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & temperature_interpolated = temperature_on_qpoints(type, ghost_type);

    // compute the temperature on quadrature points
    // this->getFEEngine().interpolateOnIntegrationPoints(
    //     *temperature, temperature_interpolated, 1, type, ghost_type);

    auto & cond = conductivity_on_qpoints(type, ghost_type);
    auto & initial_cond = initial_conductivity_array(type, ghost_type);
    auto & T_ref = T_ref_array(type, ghost_type);
    auto & cond_var = conductivity_variation_array(type, ghost_type);
    for (auto && tuple :
         zip(make_view(cond, spatial_dimension, spatial_dimension),
             temperature_interpolated,
             make_view(initial_cond, spatial_dimension, spatial_dimension),
             T_ref, cond_var)) {
      auto & C = std::get<0>(tuple);
      auto & T = std::get<1>(tuple);
      auto & init_C = std::get<2>(tuple);
      auto & T_reference = std::get<3>(tuple);
      auto & C_var = std::get<4>(tuple);
      C = init_C;

      Matrix<Real> variation(spatial_dimension, spatial_dimension,
                             C_var * (T - T_reference));
      C += variation;
    }
  }

  conductivity_release[ghost_type] = temperature_release;
  initial_conductivity[ghost_type] = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeKgradT(GhostType ghost_type) {
  computeConductivityOnQuadPoints(ghost_type);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & gradient = temperature_gradient(type, ghost_type);
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
void HeatTransferModel::computeGradT(GhostType ghost_type) {
  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & gradient = temperature_gradient(type, ghost_type);
    this->getFEEngine().gradientOnIntegrationPoints(*temperature, gradient, 1,
                                                    type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleInternalHeatRate() {
  AKANTU_DEBUG_IN();

  this->internal_heat_rate->zero();

  computeGradT(_not_ghost);

  // communicate the stresses
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  this->asynchronousSynchronize(SynchronizationTag::_htm_gradient_temperature);

  assembleInternalHeatRate(_not_ghost);

  // finalize communications
  AKANTU_DEBUG_INFO("Wait distant stresses");
  this->waitEndSynchronize(SynchronizationTag::_htm_gradient_temperature);

  assembleInternalHeatRate(_ghost);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleInternalHeatRate(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  auto & fem = this->getFEEngine();
  computeKgradT(ghost_type);

  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    auto & k_gradt_on_qpoints_vect = k_gradt_on_qpoints(type, ghost_type);

    UInt nb_quad_points = k_gradt_on_qpoints_vect.size();
    Array<Real> bt_k_gT(nb_quad_points, nb_nodes_per_element);
    fem.computeBtD(k_gradt_on_qpoints_vect, bt_k_gT, type, ghost_type);

    UInt nb_elements = mesh.getNbElement(type, ghost_type);
    Array<Real> int_bt_k_gT(nb_elements, nb_nodes_per_element);

    fem.integrate(bt_k_gT, int_bt_k_gT, nb_nodes_per_element, type, ghost_type);

    this->getDOFManager().assembleElementalArrayLocalArray(
        int_bt_k_gT, *this->internal_heat_rate, type, ghost_type, 1);
  }
}

/* --------------------------------------------------------------------------
 */
void HeatTransferModel::computeTempOnQpoints(GhostType ghost_type) {

  auto & fem = this->getFEEngine();

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & t_on_qpoints = temperature_on_qpoints(type, ghost_type);

    // compute the temperature on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        *temperature, t_on_qpoints, 1, type, ghost_type);
  }
  this->synchronizeField(SynchronizationTag::_htm_temperature_on_qpoints);
}
/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real conductivitymax = conductivity(0, 0);
  Real densitymin = density;
  Real capacitymin = capacity;
  // get maximum conductivity, and minimum density and capacity
  for (auto type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
    for (auto && data : zip(make_view(conductivity_on_qpoints(type),
                                      spatial_dimension, spatial_dimension),
                            density_array(type), capacity_array(type))) {
      auto && conductivity_ip = std::get<0>(data);
      auto && density_ip = std::get<1>(data);
      auto && capacity_ip = std::get<2>(data);

      // get the biggest parameter from k11 until k33//
      for (UInt i = 0; i < spatial_dimension; i++) {
        for (UInt j = 0; j < spatial_dimension; j++) {
          conductivitymax = std::max(conductivity_ip(i, j), conductivitymax);
        }
      }
      // get the smallest density and capacity
      densitymin = std::min(density_ip, densitymin);
      capacitymin = std::min(capacity_ip, capacitymin);
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
                      << ", the max conductivity is : " << conductivitymax
                      << ", the min density is : " << densitymin
                      << ", and the min capacity is : " << capacitymin);
  }

  Real min_dt = min_el_size * min_el_size / 4. * densitymin * capacitymin /
                conductivitymax;

  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}
/* -------------------------------------------------------------------------- */

void HeatTransferModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

  this->mesh.getDumper("heat_transfer").setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::readMaterials() {
  auto sect = this->getParserSection();

  if (not std::get<1>(sect)) {
    const auto & section = std::get<0>(sect);
    this->parseSection(section);
  }

  conductivity_on_qpoints.set(conductivity);
  initial_conductivity_array.set(conductivity);
  density_array.set(density);
  capacity_array.set(capacity);
  conductivity_variation_array.set(conductivity_variation);
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  readMaterials();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacity() {
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
void HeatTransferModel::computeRho(Array<Real> & rho, ElementType type,
                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = this->getFEEngine();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  rho.resize(nb_element * nb_quadrature_points);
  for (auto && data : zip(rho, capacity_array(type, ghost_type),
                          density_array(type, ghost_type))) {
    auto && rho_ip = std::get<0>(data);
    auto && capacity_ip = std::get<1>(data);
    auto && density_ip = std::get<2>(data);

    rho_ip = capacity_ip * density_ip;
  }

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
Real HeatTransferModel::computeThermalEnergyByNode() {
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
        heat += heat_rate[i] * this->getTimeStep();
      }
    }
    ethermal += heat;
  }

  mesh.getCommunicator().allReduce(ethermal, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return ethermal;
}

/* -------------------------------------------------------------------------- */
template <class iterator>
void HeatTransferModel::getThermalEnergy(
    iterator Eth, Array<Real>::const_iterator<Real> T_it,
    const Array<Real>::const_iterator<Real> & T_end,
    Array<Real>::const_iterator<Real> capacity_it,
    Array<Real>::const_iterator<Real> density_it) const {
  for (; T_it != T_end; ++T_it, ++Eth, ++capacity_it, ++density_it) {
    *Eth = *capacity_it * *density_it * *T_it;
  }
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy(ElementType type, UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);
  Vector<Real> Eth_on_quarature_points(nb_quadrature_points);

  auto T_it = this->temperature_on_qpoints(type).begin();
  T_it += index * nb_quadrature_points;
  auto capacity_it = capacity_array(type, _not_ghost).begin();
  capacity_it += index * nb_quadrature_points;
  auto density_it = density_array(type, _not_ghost).begin();
  density_it += index * nb_quadrature_points;

  auto T_end = T_it + nb_quadrature_points;

  getThermalEnergy(Eth_on_quarature_points.storage(), T_it, T_end, capacity_it,
                   density_it);

  return getFEEngine().integrate(Eth_on_quarature_points, type, index);
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy() {
  Real Eth = 0;

  auto & fem = getFEEngine();

  computeTempOnQpoints(_not_ghost);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
    auto nb_element = mesh.getNbElement(type, _not_ghost);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, _not_ghost);
    Array<Real> Eth_per_quad(nb_element * nb_quadrature_points, 1);

    auto & temperature_interpolated = temperature_on_qpoints(type);

    // compute the temperature on quadrature points
    // this->getFEEngine().interpolateOnIntegrationPoints(
    //     *temperature, temperature_interpolated, 1, type);

    auto T_it = temperature_interpolated.begin();
    auto T_end = temperature_interpolated.end();
    auto capacity_it = capacity_array(type, _not_ghost).begin();
    auto density_it = density_array(type, _not_ghost).begin();
    getThermalEnergy(Eth_per_quad.begin(), T_it, T_end, capacity_it,
                     density_it);

    Eth += fem.integrate(Eth_per_quad, type);
  }

  return Eth;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id) {
  AKANTU_DEBUG_IN();
  Real energy = 0;

  if (id == "thermal") {
    energy = getThermalEnergy();
  }
  // reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id, ElementType type,
                                  UInt index) {
  AKANTU_DEBUG_IN();

  Real energy = 0.;

  if (id == "thermal") {
    energy = getThermalEnergy(type, index);
  }

  AKANTU_DEBUG_OUT();
  return energy;
}
/* -------------------------------------------------------------------------- */
void HeatTransferModel::assignPropertyToPhysicalGroup(
    const std::string & property_name, const std::string & group_name,
    Real value) {
  AKANTU_DEBUG_ASSERT(property_name != "conductivity",
                      "Scalar was provided instead of a conductivity matrix");
  auto && el_group = mesh.getElementGroup(group_name);
  auto & fem = this->getFEEngine();

  for (auto ghost_type : ghost_types) {
    for (auto && type :
         mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
      auto nb_quadrature_points = fem.getNbIntegrationPoints(type);
      auto && elements = el_group.getElements(type, ghost_type);

      Array<Real>::scalar_iterator field_it;
      if (property_name == "density") {
        field_it = density_array(type, ghost_type).begin();
      } else if (property_name == "capacity") {
        field_it = capacity_array(type, ghost_type).begin();
      } else {
        AKANTU_EXCEPTION(property_name +
                         " is not a valid material property name.");
      }
      for (auto && el : elements) {
        for (auto && qpoint : arange(nb_quadrature_points)) {
          field_it[el * nb_quadrature_points + qpoint] = value;
        }
      }
      need_to_reassemble_capacity = true;
      need_to_reassemble_capacity_lumped = true;
    }
  }
}
/* -------------------------------------------------------------------------- */
void HeatTransferModel::assignPropertyToPhysicalGroup(
    const std::string & property_name, const std::string & group_name,
    Matrix<Real> cond_matrix) {
  AKANTU_DEBUG_ASSERT(property_name == "conductivity",
                      "When updating material parameters, only conductivity "
                      "accepts matrix as an input");
  auto && el_group = mesh.getElementGroup(group_name);
  auto & fem = this->getFEEngine();
  auto dim = this->getMesh().getSpatialDimension();

  for (auto ghost_type : ghost_types) {
    for (auto && type :
         mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
      auto nb_quadrature_points = fem.getNbIntegrationPoints(type);
      auto && elements = el_group.getElements(type, ghost_type);

      auto init_cond_it =
          make_view(initial_conductivity_array(type, ghost_type), dim, dim)
              .begin();
      auto cond_on_quad_it =
          make_view(conductivity_on_qpoints(type, ghost_type), dim, dim)
              .begin();

      for (auto && el : elements) {
        for (auto && qpoint : arange(nb_quadrature_points)) {
          init_cond_it[el * nb_quadrature_points + qpoint] = cond_matrix;
          cond_on_quad_it[el * nb_quadrature_points + qpoint] = cond_matrix;
        }
      }
    }
    conductivity_release[ghost_type] += 1;
  }
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> HeatTransferModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  auto field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> HeatTransferModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION(
        "Capacity lumped is a nodal field now stored in the DOF manager."
        "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["temperature"] = temperature.get();
  real_nodal_fields["temperature_rate"] = temperature_rate.get();
  real_nodal_fields["external_heat_rate"] = external_heat_rate.get();
  real_nodal_fields["internal_heat_rate"] = internal_heat_rate.get();
  real_nodal_fields["increment"] = increment.get();

  std::shared_ptr<dumpers::Field> field =
      mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> HeatTransferModel::createElementalField(
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
  } else if (field_name == "temperature_on_qpoints") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(temperature_on_qpoints);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        temperature_on_qpoints, group_name, this->spatial_dimension,
        element_kind, nb_data_per_elem);
  } else if (field_name == "conductivity") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(conductivity_on_qpoints);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        conductivity_on_qpoints, group_name, this->spatial_dimension,
        element_kind, nb_data_per_elem);
  } else if (field_name == "density") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(density_array);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        density_array, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  } else if (field_name == "capacity") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(capacity_array);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        capacity_array, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
