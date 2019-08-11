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
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of HeatTransferModel class
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
#include "fluid_diffusion_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
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

namespace fluid_diffusion {
  namespace details {
    class ComputeRhoFunctor {
    public:
      ComputeRhoFunctor(const FluidDiffusionModel & model) : model(model){};

      void operator()(Matrix<Real> & rho, const Element &) const {
        rho.set(model.getCapacity() * model.getDensity());
      }

    private:
      const FluidDiffusionModel & model;
    };
  } // namespace details
} // namespace fluid_diffusion

/* -------------------------------------------------------------------------- */
FluidDiffusionModel::FluidDiffusionModel(Mesh & mesh, UInt dim, const ID & id,
                                         const MemoryID & memory_id)
    : Model(mesh, ModelType::_fluid_diffusion_model, dim, id, memory_id),
      pressure_gradient("pressure_gradient", id),
      pressure_on_qpoints("pressure_on_qpoints", id),
      permeability_on_qpoints("permeability_on_qpoints", id),
      k_gradp_on_qpoints("k_gradp_on_qpoints", id) {
  AKANTU_DEBUG_IN();

  this->initDOFManager();

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_htm_pressure);
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_gradient_pressure);
  }

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh,
                                       spatial_dimension - 1);

#ifdef AKANTU_USE_IOHELPER
  this->mesh.registerDumper<DumperParaview>("fluid_diffusion", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension - 1, _not_ghost, _ek_regular);
#endif

  this->registerParam("viscosity", viscosity, _pat_parsmod);
  this->registerParam("pressure_reference", P_ref, 0., _pat_parsmod);
  this->registerParam("capacity", capacity, _pat_parsmod);
  this->registerParam("density", density, _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  pressure_on_qpoints.initialize(fem, _nb_component = 1);
  pressure_gradient.initialize(fem, _nb_component = 1);
  permeability_on_qpoints.initialize(fem, _nb_component = 1);
  k_gradp_on_qpoints.initialize(fem, _nb_component = 1);
}

/* -------------------------------------------------------------------------- */
FEEngine & FluidDiffusionModel::getFEEngineBoundary(const ID & name) {
  return aka::as_type<FEEngine>(getFEEngineClassBoundary<FEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
template <typename T>
void FluidDiffusionModel::allocNodalField(Array<T> *& array, const ID & name) {
  if (array == nullptr) {
    UInt nb_nodes = mesh.getNbNodes();
    std::stringstream sstr_disp;
    sstr_disp << id << ":" << name;

    array = &(alloc<T>(sstr_disp.str(), nb_nodes, 1, T()));
  }
}

/* -------------------------------------------------------------------------- */
FluidDiffusionModel::~FluidDiffusionModel() = default;

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleCapacityLumped(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<FEEngineType>();
  fluid_diffusion::details::ComputeRhoFunctor compute_rho(*this);

  for (auto & type :
       mesh.elementTypes(spatial_dimension - 1, ghost_type, _ek_regular)) {
    fem.assembleFieldLumped(compute_rho, "M", "pressure", this->getDOFManager(),
                            type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MatrixType FluidDiffusionModel::getMatrixType(const ID & matrix_id) {
  if (matrix_id == "K" or matrix_id == "M") {
    return _symmetric;
  }

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assemblePermeabilityMatrix();
  } else if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacity();
  }
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacityLumped();
  }
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  this->assembleInternalFlux();

  this->getDOFManager().assembleToResidual("pressure", *this->external_flux, 1);
  this->getDOFManager().assembleToResidual("pressure", *this->internal_flux, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::predictor() { ++pressure_release; }

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().clearLumpedMatrix("M");

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);

  need_to_reassemble_capacity_lumped = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::initSolver(TimeStepSolverType time_step_solver_type,
                                     NonLinearSolverType) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->pressure, "pressure");
  this->allocNodalField(this->external_flux, "external_flux");
  this->allocNodalField(this->internal_flux, "internal_flux");
  this->allocNodalField(this->blocked_dofs, "blocked_dofs");

  if (!dof_manager.hasDOFs("pressure")) {
    dof_manager.registerDOFs("pressure", *this->pressure, _dst_nodal);
    dof_manager.registerBlockedDOFs("pressure", *this->blocked_dofs);
  }

  // if (time_step_solver_type == TimeStepSolverType::_dynamic ||
  //     time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
  //   this->allocNodalField(this->pressure_rate, "pressure_rate");

  if (!dof_manager.hasDOFsDerivatives("pressure", 1)) {
    dof_manager.registerDOFsDerivative("pressure", 1, *this->pressure_rate);
  }
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
FluidDiffusionModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  // case _explicit_lumped_mass: {
  //   return std::make_tuple("explicit_lumped",
  //   TimeStepSolverType::_dynamic_lumped);
  // }
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  // case _implicit_dynamic: {
  //   return std::make_tuple("implicit", TimeStepSolverType::_dynamic);
  // }
  default:
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }
}

/* --------------------------------------------------------------------------
 */
ModelSolverOptions FluidDiffusionModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  // case TimeStepSolverType::_dynamic_lumped: {
  //   options.non_linear_solver_type = NonLinearSolverType::_lumped;
  //   options.integration_scheme_type["pressure"] =
  //   IntegrationSchemeType::_forward_euler;
  //   options.solution_type["pressure"] = IntegrationScheme::_pressure_rate;
  //   break;
  // }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["pressure"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["pressure"] = IntegrationScheme::_not_defined;
    break;
  }
  // case TimeStepSolverType::_dynamic: {
  //   if (this->method == _explicit_consistent_mass) {
  //     options.non_linear_solver_type =
  //     NonLinearSolverType::_newton_raphson;
  //     options.integration_scheme_type["pressure"] =
  //     IntegrationSchemeType::_forward_euler;
  //     options.solution_type["pressure"] =
  //         IntegrationScheme::_pressure_rate;
  //   } else {
  //     options.non_linear_solver_type =
  //     NonLinearSolverType::_newton_raphson;
  //     options.integration_scheme_type["pressure"] =
  //     IntegrationSchemeType::_backward_euler;
  //     options.solution_type["pressure"] = IntegrationScheme::_pressure;
  //   }
  //   break;
  // }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* --------------------------------------------------------------------------
 */
void FluidDiffusionModel::assemblePermeabilityMatrix() {
  AKANTU_DEBUG_IN();

  this->computePermeabilityOnQuadPoints(_not_ghost);
  if (permeability_release[_not_ghost] == permeability_matrix_release)
    return;

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }
  this->getDOFManager().clearMatrix("K");

  this->assemblePermeabilityMatrix(_not_ghost);
  // switch (mesh.getSpatialDimension()) {
  // case 1:
  // this->assemblePermeabilityMatrix<1>(_not_ghost);
  //   break;
  // case 2:
  // this->assemblePermeabilityMatrix<2>(_not_ghost);
  //   break;
  // case 3:
  //   this->assemblePermeabilityMatrix<3>(_not_ghost);
  //   break;
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assemblePermeabilityMatrix(
    const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();

  for (auto && type :
       mesh.elementTypes(spatial_dimension - 1, ghost_type, _ek_regular)) {
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto bt_d_b = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    fem.computeBtDB(permeability_on_qpoints(type, ghost_type), *bt_d_b, 2, type,
                    ghost_type);

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    auto K_e = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_e");

    fem.integrate(*bt_d_b, *K_e, nb_nodes_per_element * nb_nodes_per_element,
                  type, ghost_type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "pressure", *K_e, type, ghost_type, _symmetric);
  }

  permeability_matrix_release = permeability_release[ghost_type];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::computePermeabilityOnQuadPoints(
    const GhostType & ghost_type) {
  // // if already computed once check if need to compute
  // if (not initial_permeability[ghost_type]) {
  //   // if pressure did not change, conductivity will not vary
  //   if (pressure_release == permeability_release[ghost_type])
  //     return;

  //   // if permeability_variation is 0 no need to recompute
  //   if (permeability_variation == 0.)
  //     return;
  // }

  for (auto & type :
       mesh.elementTypes(spatial_dimension - 1, ghost_type, _ek_regular)) {

    auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);
    auto nb_element = mesh.getNbElement(type);
    Array<Real> apperture_interpolated(nb_element * nb_nodes_per_element, 1);

    // compute the pressure on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        *apperture, apperture_interpolated, 1, type, ghost_type);

    auto & perm = permeability_on_qpoints(type, ghost_type);
    for (auto && tuple : zip(perm, apperture_interpolated)) {
      auto & k = std::get<0>(tuple);
      auto & apper = std::get<1>(tuple);
      k = apper * apper / (12 * this->viscosity);
    }
  }

  permeability_release[ghost_type] = pressure_release;
  initial_permeability[ghost_type] = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::computeKgradP(const GhostType & ghost_type) {
  computePermeabilityOnQuadPoints(ghost_type);

  for (auto & type :
       mesh.elementTypes(spatial_dimension - 1, ghost_type, _ek_regular)) {
    auto & gradient = pressure_gradient(type, ghost_type);
    this->getFEEngine().gradientOnIntegrationPoints(*pressure, gradient, 1,
                                                    type, ghost_type);

    for (auto && values : zip(permeability_on_qpoints(type, ghost_type),
                              gradient, k_gradp_on_qpoints(type, ghost_type))) {
      const auto & k = std::get<0>(values);
      const auto & BP = std::get<1>(values);
      auto & k_BP = std::get<2>(values);

      k_BP = k * BP;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleInternalFlux() {
  AKANTU_DEBUG_IN();

  this->internal_flux->clear();

  this->synchronize(SynchronizationTag::_htm_pressure);
  auto & fem = this->getFEEngine();

  for (auto ghost_type : ghost_types) {
    // compute k \grad P
    computeKgradP(ghost_type);

    for (auto type :
         mesh.elementTypes(spatial_dimension - 1, ghost_type, _ek_regular)) {
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      auto & k_gradp_on_qpoints_vect = k_gradp_on_qpoints(type, ghost_type);

      UInt nb_quad_points = k_gradp_on_qpoints_vect.size();
      Array<Real> bt_k_gP(nb_quad_points, nb_nodes_per_element);
      fem.computeBtD(k_gradp_on_qpoints_vect, bt_k_gP, type, ghost_type);

      UInt nb_elements = mesh.getNbElement(type, ghost_type);
      Array<Real> int_bt_k_gP(nb_elements, nb_nodes_per_element);

      fem.integrate(bt_k_gP, int_bt_k_gP, nb_nodes_per_element, type,
                    ghost_type);

      this->getDOFManager().assembleElementalArrayLocalArray(
          int_bt_k_gP, *this->internal_flux, type, ghost_type, -1);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
// Real FluidDiffusionModel::getStableTimeStep() {
//   AKANTU_DEBUG_IN();

//   Real el_size;
//   Real min_el_size = std::numeric_limits<Real>::max();
//   Real permeabilitymax = 1 / this->viscosity;

//   for (auto & type :
//        mesh.elementTypes(spatial_dimension - 1, _not_ghost, _ek_regular)) {

//     UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

//     Array<Real> coord(0, nb_nodes_per_element * spatial_dimension);
//     FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), coord, type,
//                                          _not_ghost);

//     auto el_coord = coord.begin(spatial_dimension, nb_nodes_per_element);
//     UInt nb_element = mesh.getNbElement(type);

//     for (UInt el = 0; el < nb_element; ++el, ++el_coord) {
//       el_size = getFEEngine().getElementInradius(*el_coord, type);
//       min_el_size = std::min(min_el_size, el_size);
//     }

//     AKANTU_DEBUG_INFO("The minimum element size : "
//                       << min_el_size
//                       << " and the max permeability is : " <<
//                       permeabilitymax);
//   }

//   Real min_dt = 2. * min_el_size * min_el_size / 4. * density * capacity /
//                 permeabilitymax;

//   mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

//   AKANTU_DEBUG_OUT();

//   return min_dt;
// }
// /* --------------------------------------------------------------------------
// */

void FluidDiffusionModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper("fluid_diffusion").setTimeStep(time_step);
#endif
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::readMaterials() {
  auto sect = this->getParserSection();

  if (not std::get<1>(sect)) {
    const auto & section = std::get<0>(sect);
    this->parseSection(section);
  }
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  readMaterials();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleCapacity() {
  AKANTU_DEBUG_IN();
  auto ghost_type = _not_ghost;

  this->getDOFManager().clearMatrix("M");

  auto & fem = getFEEngineClass<FEEngineType>();

  fluid_diffusion::details::ComputeRhoFunctor rho_functor(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension - 1, ghost_type, _ek_regular)) {
    fem.assembleFieldMatrix(rho_functor, "M", "pressure", this->getDOFManager(),
                            type, ghost_type);
  }

  need_to_reassemble_capacity = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::computeRho(Array<Real> & rho, ElementType type,
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

// /* --------------------------------------------------------------------------
// */ Real FluidDiffusionModel::computeThermalEnergyByNode() {
//   AKANTU_DEBUG_IN();

//   Real ethermal = 0.;

//   for (auto && pair : enumerate(make_view(
//            *internal_heat_rate, internal_heat_rate->getNbComponent()))) {
//     auto n = std::get<0>(pair);
//     auto & heat_rate = std::get<1>(pair);

//     Real heat = 0.;
//     bool is_local_node = mesh.isLocalOrMasterNode(n);
//     bool count_node = is_local_node;

//     for (UInt i = 0; i < heat_rate.size(); ++i) {
//       if (count_node)
//         heat += heat_rate[i] * time_step;
//     }
//     ethermal += heat;
//   }

//   mesh.getCommunicator().allReduce(ethermal, SynchronizerOperation::_sum);

//   AKANTU_DEBUG_OUT();
//   return ethermal;
// }

// /* --------------------------------------------------------------------------
// */ template <class iterator> void FluidDiffusionModel::getThermalEnergy(
//     iterator Eth, Array<Real>::const_iterator<Real> T_it,
//     Array<Real>::const_iterator<Real> T_end) const {
//   for (; T_it != T_end; ++T_it, ++Eth) {
//     *Eth = capacity * density * *T_it;
//   }
// }

// /* --------------------------------------------------------------------------
// */ Real FluidDiffusionModel::getThermalEnergy(const ElementType & type, UInt
// index) {
//   AKANTU_DEBUG_IN();

//   UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);
//   Vector<Real> Eth_on_quarature_points(nb_quadrature_points);

//   auto T_it = this->pressure_on_qpoints(type).begin();
//   T_it += index * nb_quadrature_points;

//   auto T_end = T_it + nb_quadrature_points;

//   getThermalEnergy(Eth_on_quarature_points.storage(), T_it, T_end);

//   return getFEEngine().integrate(Eth_on_quarature_points, type, index);
// }

// /* --------------------------------------------------------------------------
// */ Real FluidDiffusionModel::getThermalEnergy() {
//   Real Eth = 0;

//   auto & fem = getFEEngine();

//   for (auto && type :
//        mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
//     auto nb_element = mesh.getNbElement(type, _not_ghost);
//     auto nb_quadrature_points = fem.getNbIntegrationPoints(type, _not_ghost);
//     Array<Real> Eth_per_quad(nb_element * nb_quadrature_points, 1);

//     auto & pressure_interpolated = pressure_on_qpoints(type);

//     // compute the pressure on quadrature points
//     this->getFEEngine().interpolateOnIntegrationPoints(
//         *pressure, pressure_interpolated, 1, type);

//     auto T_it = pressure_interpolated.begin();
//     auto T_end = pressure_interpolated.end();
//     getThermalEnergy(Eth_per_quad.begin(), T_it, T_end);

//     Eth += fem.integrate(Eth_per_quad, type);
//   }

//   return Eth;
// }

// /* --------------------------------------------------------------------------
// */ Real FluidDiffusionModel::getEnergy(const std::string & id) {
//   AKANTU_DEBUG_IN();
//   Real energy = 0;

//   if (id == "thermal")
//     energy = getThermalEnergy();

//   // reduction sum over all processors
//   mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

//   AKANTU_DEBUG_OUT();
//   return energy;
// }

// /* --------------------------------------------------------------------------
// */ Real FluidDiffusionModel::getEnergy(const std::string & id,
//                                   const ElementType & type, UInt index) {
//   AKANTU_DEBUG_IN();

//   Real energy = 0.;

//   if (id == "thermal")
//     energy = getThermalEnergy(type, index);

//   AKANTU_DEBUG_OUT();
//   return energy;
// }

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

std::shared_ptr<dumper::Field> FluidDiffusionModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs;

  auto field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field> FluidDiffusionModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION(
        "Capacity lumped is a nodal field now stored in the DOF manager."
        "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["pressure"] = pressure;
  real_nodal_fields["pressure_rate"] = pressure_rate;
  real_nodal_fields["external_flux"] = external_flux;
  real_nodal_fields["internal_flux"] = internal_flux;
  real_nodal_fields["increment"] = increment;
  real_nodal_fields["apperture"] = apperture;

  std::shared_ptr<dumper::Field> field =
      mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field> FluidDiffusionModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const UInt & spatial_dimension,
    const ElementKind & element_kind) {

  std::shared_ptr<dumper::Field> field;

  if (field_name == "partitions")
    field = mesh.createElementalField<UInt, dumper::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension - 1,
        element_kind);
  else if (field_name == "pressure_gradient") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(pressure_gradient, element_kind);

    field = mesh.createElementalField<Real, dumper::InternalMaterialField>(
        pressure_gradient, group_name, this->spatial_dimension - 1,
        element_kind, nb_data_per_elem);
  } else if (field_name == "permeability") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(permeability_on_qpoints, element_kind);

    field = mesh.createElementalField<Real, dumper::InternalMaterialField>(
        permeability_on_qpoints, group_name, this->spatial_dimension - 1,
        element_kind, nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field> FluidDiffusionModel::createElementalField(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const ElementKind & element_kind) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field> FluidDiffusionModel::createNodalFieldBool(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field> FluidDiffusionModel::createNodalFieldReal(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}
#endif

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::dump(const std::string & dumper_name) {
  mesh.dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::dump(const std::string & dumper_name, UInt step) {
  mesh.dump(dumper_name, step);
}

/* ------------------------------------------------------------------------- */
void FluidDiffusionModel::dump(const std::string & dumper_name, Real time,
                               UInt step) {
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::dump() { mesh.dump(); }

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::dump(UInt step) { mesh.dump(step); }

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::dump(Real time, UInt step) { mesh.dump(time, step); }

/* -------------------------------------------------------------------------- */
inline UInt
FluidDiffusionModel::getNbData(const Array<UInt> & indexes,
                               const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = indexes.size();

  switch (tag) {
  case SynchronizationTag::_htm_pressure: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void
FluidDiffusionModel::packData(CommunicationBuffer & buffer,
                              const Array<UInt> & indexes,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_htm_pressure: {
      buffer << (*pressure)(index);
      break;
    }
    default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void FluidDiffusionModel::unpackData(CommunicationBuffer & buffer,
                                            const Array<UInt> & indexes,
                                            const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_htm_pressure: {
      buffer >> (*pressure)(index);
      break;
    }
    default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt
FluidDiffusionModel::getNbData(const Array<Element> & elements,
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
  case SynchronizationTag::_htm_pressure: {
    size += nb_nodes_per_element * sizeof(Real); // pressure
    break;
  }
  case SynchronizationTag::_htm_gradient_pressure: {
    // pressure gradient
    size += getNbIntegrationPoints(elements) * sizeof(Real);
    size += nb_nodes_per_element * sizeof(Real); // nodal pressures
    break;
  }
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void
FluidDiffusionModel::packData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_htm_pressure: {
    packNodalDataHelper(*pressure, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_htm_gradient_pressure: {
    packElementalDataHelper(pressure_gradient, buffer, elements, true,
                            getFEEngine());
    packNodalDataHelper(*pressure, buffer, elements, mesh);
    break;
  }
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }
}

/* -------------------------------------------------------------------------- */
inline void FluidDiffusionModel::unpackData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_htm_pressure: {
    unpackNodalDataHelper(*pressure, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_htm_gradient_pressure: {
    unpackElementalDataHelper(pressure_gradient, buffer, elements, true,
                              getFEEngine());
    unpackNodalDataHelper(*pressure, buffer, elements, mesh);

    break;
  }
  default: { AKANTU_ERROR("Unknown ghost synchronization tag : " << tag); }
  }
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
