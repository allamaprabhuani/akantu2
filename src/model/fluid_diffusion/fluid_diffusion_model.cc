/**
 * @file   fluid_diffusion_model.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Aug 2019
 *
 * @brief  Implementation of Transient fluid diffusion model
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
#include "group_manager_inline_impl.hh"
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

      void operator()(Matrix<Real> & rho, const Element & element) const {
        const auto & type = element.type;
        const auto & gt = element.ghost_type;
        const auto & aperture = model.getApertureOnQpoints(type, gt);
        auto nb_quad_points = aperture.getNbComponent();
        Real aperture_av = 0.;
        for (auto q : arange(nb_quad_points)) {
          aperture_av += aperture(element.element, q);
        }
        aperture_av /= nb_quad_points;
        if (model.isModelVelocityDependent()) {
          rho.set(model.getCompressibility() * aperture_av);
        } else {
          rho.set((model.getCompressibility() + model.getPushability()) *
                  aperture_av);
        }
      }

    private:
      const FluidDiffusionModel & model;
    };
  } // namespace details
} // namespace fluid_diffusion

/* -------------------------------------------------------------------------- */
FluidDiffusionModel::FluidDiffusionModel(Mesh & mesh, UInt dim, const ID & id)
    : Model(mesh, ModelType::_fluid_diffusion_model, dim, id),
      pressure_gradient("pressure_gradient", id),
      pressure_on_qpoints("pressure_on_qpoints", id),
      delta_pres_on_qpoints("delta_pres_on_qpoints", id),
      permeability_on_qpoints("permeability_on_qpoints", id),
      aperture_on_qpoints("aperture_on_qpoints", id),
      prev_aperture_on_qpoints("prev_aperture_on_qpoints", id),
      k_gradp_on_qpoints("k_gradp_on_qpoints", id) {
  AKANTU_DEBUG_IN();

  // mesh.registerEventHandler(*this, akantu::_ehp_lowest);

  this->spatial_dimension = std::max(1, int(mesh.getSpatialDimension()) - 1);

  this->initDOFManager();

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_fdm_pressure);
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_fdm_gradient_pressure);
  }

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh, spatial_dimension);

#ifdef AKANTU_USE_IOHELPER
  this->mesh.registerDumper<DumperParaview>("fluid_diffusion", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif

  this->registerParam("viscosity", viscosity, 1., _pat_parsmod);
  this->registerParam("compressibility", compressibility, 1., _pat_parsmod);
  this->registerParam("default_aperture", default_aperture, 1.e-5,
                      _pat_parsmod);
  this->registerParam("use_aperture_speed", use_aperture_speed, _pat_parsmod);
  this->registerParam("pushability", pushability, 1., _pat_parsmod);
  this->registerParam("insertion_damage", insertion_damage, 0., _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  pressure_on_qpoints.initialize(fem, _nb_component = 1);
  delta_pres_on_qpoints.initialize(fem, _nb_component = 1);
  pressure_gradient.initialize(fem, _nb_component = 1);
  permeability_on_qpoints.initialize(fem, _nb_component = 1);
  aperture_on_qpoints.initialize(fem, _nb_component = 1);
  if (use_aperture_speed) {
    prev_aperture_on_qpoints.initialize(fem, _nb_component = 1);
  }
  k_gradp_on_qpoints.initialize(fem, _nb_component = 1);
}

/* -------------------------------------------------------------------------- */
FEEngine & FluidDiffusionModel::getFEEngineBoundary(const ID & name) {
  return aka::as_type<FEEngine>(getFEEngineClassBoundary<FEEngineType>(name));
}
/* -------------------------------------------------------------------------- */
FluidDiffusionModel::~FluidDiffusionModel() = default;

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleCapacityLumped(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<FEEngineType>();
  fluid_diffusion::details::ComputeRhoFunctor compute_rho(*this);

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
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
void FluidDiffusionModel::afterSolveStep(bool converged) {
  AKANTU_DEBUG_IN();

  if (not converged) {
    return;
  }

  ++solution_release;
  need_to_reassemble_capacity_lumped = true;
  need_to_reassemble_capacity = true;

  /// interpolating pressures on integration points for further use by BC
  const GhostType gt = _not_ghost;
  for (auto && type : mesh.elementTypes(spatial_dimension, gt, _ek_regular)) {
    Array<Real> tmp(pressure_on_qpoints(type, gt));
    this->getFEEngine().interpolateOnIntegrationPoints(
        *pressure, pressure_on_qpoints(type, gt), 1, type, gt);
    delta_pres_on_qpoints(type, gt) = pressure_on_qpoints(type, gt);
    delta_pres_on_qpoints(type, gt) -= tmp;
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().getLumpedMatrix("M").zero();

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);

  need_to_reassemble_capacity_lumped = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::initSolver(
    TimeStepSolverType time_step_solver_type,
    NonLinearSolverType /*non_linear_solver_type*/) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->pressure, 1, "pressure");
  this->allocNodalField(this->external_flux, 1, "external_flux");
  this->allocNodalField(this->internal_flux, 1, "internal_flux");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");
  // this->allocNodalField(this->aperture, "aperture");

  if (!dof_manager.hasDOFs("pressure")) {
    dof_manager.registerDOFs("pressure", *this->pressure, _dst_nodal);
    dof_manager.registerBlockedDOFs("pressure", *this->blocked_dofs);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->pressure_rate, 1, "pressure_rate");

    if (!dof_manager.hasDOFsDerivatives("pressure", 1)) {
      dof_manager.registerDOFsDerivative("pressure", 1, *this->pressure_rate);
    }
  }
}
/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
FluidDiffusionModel::getDefaultSolverID(const AnalysisMethod & method) {
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
ModelSolverOptions FluidDiffusionModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["pressure"] =
        IntegrationSchemeType::_forward_euler;
    options.solution_type["pressure"] = IntegrationScheme::_pressure_rate;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["pressure"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["pressure"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["pressure"] =
          IntegrationSchemeType::_forward_euler;
      options.solution_type["pressure"] = IntegrationScheme::_pressure_rate;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["pressure"] =
          IntegrationSchemeType::_backward_euler;
      options.solution_type["pressure"] = IntegrationScheme::_pressure;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assemblePermeabilityMatrix() {
  AKANTU_DEBUG_IN();

  this->computePermeabilityOnQuadPoints(_not_ghost);
  if (permeability_release[_not_ghost] == permeability_matrix_release) {
    return;
  }

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }
  this->getDOFManager().getMatrix("K").zero();

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
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
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

  // if already computed once check if need to compute
  if (not initial_permeability[ghost_type]) {
    // if aperture did not change, permeability will not vary
    if (solution_release == permeability_release[ghost_type]) {
      return;
    }
  }

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {

    // auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);
    // auto nb_element = mesh.getNbElement(type);
    // Array<Real> aperture_interpolated(nb_element * nb_nodes_per_element,
    // 1);

    // // compute the pressure on quadrature points
    // this->getFEEngine().interpolateOnIntegrationPoints(
    //     *aperture, aperture_interpolated, 1, type, ghost_type);

    auto & perm = permeability_on_qpoints(type, ghost_type);
    auto & aperture = aperture_on_qpoints(type, ghost_type);
    for (auto && tuple : zip(perm, aperture)) {
      auto & k = std::get<0>(tuple);
      auto & aper = std::get<1>(tuple);
      k = aper * aper * aper / (12 * this->viscosity);
    }
  }

  permeability_release[ghost_type] = solution_release;
  initial_permeability[ghost_type] = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::computeKgradP(const GhostType & ghost_type) {
  computePermeabilityOnQuadPoints(ghost_type);

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
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

  this->synchronize(SynchronizationTag::_fdm_pressure);
  auto & fem = this->getFEEngine();

  for (auto ghost_type : ghost_types) {
    // compute k \grad P
    computeKgradP(ghost_type);

    for (auto type :
         mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      auto & k_gradp_on_qpoints_vect = k_gradp_on_qpoints(type, ghost_type);

      UInt nb_quad_points = k_gradp_on_qpoints_vect.size();
      Array<Real> bt_k_gP(nb_quad_points, nb_nodes_per_element);
      fem.computeBtD(k_gradp_on_qpoints_vect, bt_k_gP, type, ghost_type);

      if (this->use_aperture_speed) {
        auto aper_speed = aperture_on_qpoints(type, ghost_type);
        aper_speed -= prev_aperture_on_qpoints(type, ghost_type);
        aper_speed *= 1 / this->time_step;
        Array<Real> aper_speed_by_shapes(nb_quad_points, nb_nodes_per_element);

        fem.computeNtb(aper_speed, aper_speed_by_shapes, type, ghost_type);
        bt_k_gP += aper_speed_by_shapes;
      }
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
Real FluidDiffusionModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real max_aper_on_qpoint = std::numeric_limits<Real>::min();

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {

    auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);
    auto & aper = aperture_on_qpoints(type, _not_ghost);
    auto nb_quad_per_element = aper.getNbComponent();
    auto mesh_dim = this->mesh.getSpatialDimension();
    Array<Real> coord(0, nb_nodes_per_element * mesh_dim);
    FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), coord, type,
                                         _not_ghost);

    for (auto && data : zip(make_view(coord, mesh_dim, nb_nodes_per_element),
                            make_view(aper, nb_quad_per_element))) {
      Matrix<Real> & el_coord = std::get<0>(data);
      Vector<Real> & el_aper = std::get<1>(data);
      el_size = getFEEngine().getElementInradius(el_coord, type);
      min_el_size = std::min(min_el_size, el_size);
      max_aper_on_qpoint = std::max(max_aper_on_qpoint, el_aper.norm<L_inf>());
    }

    AKANTU_DEBUG_INFO("The minimum element size : "
                      << min_el_size
                      << " and the max aperture is : " << max_aper_on_qpoint);
  }

  Real min_dt = 2. * min_el_size * min_el_size / 4 * this->compressibility /
                max_aper_on_qpoint / max_aper_on_qpoint;

  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::setTimeStep(Real dt, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);
  this->time_step = dt;
#if defined(AKANTU_USE_IOHELPER)
  // this->mesh.getDumper("fluid_diffusion").setTimeStep(time_step);
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
  readMaterials();
  Model::initFullImpl(options);
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::assembleCapacity() {
  AKANTU_DEBUG_IN();
  auto ghost_type = _not_ghost;

  this->getDOFManager().getMatrix("M").zero();

  auto & fem = getFEEngineClass<FEEngineType>();

  fluid_diffusion::details::ComputeRhoFunctor rho_functor(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
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

  rho.copy(aperture_on_qpoints(type, ghost_type));
  if (this->use_aperture_speed) {
    rho *= this->compressibility;
  } else {
    rho *= (this->compressibility + this->pushability);
  }

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::onElementsAdded(const Array<Element> & element_list,
                                          const NewElementsEvent & /*event*/) {
  AKANTU_DEBUG_IN();
  if (element_list.empty()) {
    return;
  }
  /// TODO had to do this ugly way because the meeting is in 3 days
  // removing it will make the assemble mass matrix fail. jacobians are not
  // updated
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  this->resizeFields();

  auto type_fluid =
      *mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular).begin();
  const auto nb_quad_elem =
      getFEEngine().getNbIntegrationPoints(type_fluid, _not_ghost);
  auto aperture_it =
      make_view(getApertureOnQpoints(type_fluid), nb_quad_elem).begin();

  /// assign default aperture to newly added nodes
  for (auto && element : element_list) {
    auto elem_id = element.element;
    for (auto ip : arange(nb_quad_elem)) {
      aperture_it[elem_id](ip) = this->default_aperture;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::resizeFields() {

  AKANTU_DEBUG_IN();

  /// elemental fields
  auto & fem = this->getFEEngine();
  pressure_on_qpoints.initialize(fem, _nb_component = 1);
  delta_pres_on_qpoints.initialize(fem, _nb_component = 1);
  pressure_gradient.initialize(fem, _nb_component = 1);
  permeability_on_qpoints.initialize(fem, _nb_component = 1);
  aperture_on_qpoints.initialize(fem, _nb_component = 1);
  if (use_aperture_speed) {
    prev_aperture_on_qpoints.initialize(fem, _nb_component = 1);
  }
  k_gradp_on_qpoints.initialize(fem, _nb_component = 1);

  /// nodal fields
  UInt nb_nodes = mesh.getNbNodes();

  if (pressure != nullptr) {
    pressure->resize(nb_nodes, 0.);
    ++solution_release;
  }
  if (external_flux != nullptr) {
    external_flux->resize(nb_nodes, 0.);
  }
  if (internal_flux != nullptr) {
    internal_flux->resize(nb_nodes, 0.);
  }
  if (blocked_dofs != nullptr) {
    blocked_dofs->resize(nb_nodes, false);
  }
  if (pressure_rate != nullptr) {
    pressure_rate->resize(nb_nodes, 0.);
  }
  need_to_reassemble_capacity_lumped = true;
  need_to_reassemble_capacity = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_COHESIVE_ELEMENT)
void FluidDiffusionModel::getApertureOnQpointsFromCohesive(
    const SolidMechanicsModelCohesive & coh_model, bool first_time) {
  AKANTU_DEBUG_IN();

  // get element types
  auto & coh_mesh = coh_model.getMesh();
  const auto dim = coh_mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto type = *coh_mesh.elementTypes(dim, gt, _ek_regular).begin();
  const ElementType type_facets = Mesh::getFacetType(type);
  const ElementType typecoh = FEEngine::getCohesiveElementType(type_facets);
  // const auto nb_coh_elem = coh_mesh.getNbElement(typecoh);
  const auto nb_quad_coh_elem = coh_model.getFEEngine("CohesiveFEEngine")
                                    .getNbIntegrationPoints(typecoh, gt);

  auto type_fluid =
      *mesh.elementTypes(spatial_dimension, gt, _ek_regular).begin();
  const auto nb_fluid_elem = mesh.getNbElement(type_fluid);
  const auto nb_quad_elem =
      getFEEngine().getNbIntegrationPoints(type_fluid, gt);
  // AKANTU_DEBUG_ASSERT(nb_quad_coh_elem == nb_quad_elem,
  //                     "Different number of integration points per cohesive "
  //                     "and fluid elements");

  // save previous aperture
  if (not first_time and use_aperture_speed) {
    prev_aperture_on_qpoints.copy(aperture_on_qpoints);
  }

  // AKANTU_DEBUG_ASSERT(nb_coh_elem == nb_fluid_elem,
  //                     "Different number of cohesive and fluid elements");

  auto aperture_it =
      make_view(getApertureOnQpoints(type_fluid), nb_quad_elem).begin();

  /// loop on each segment element
  for (auto && element_nb : arange(nb_fluid_elem)) {
    Element element{type_fluid, element_nb, gt};
    auto global_coh_elem =
        mesh.getElementalData<akantu::Element>("cohesive_elements")(element);
    auto ge = global_coh_elem.element;
    auto le = coh_model.getMaterialLocalNumbering(
        global_coh_elem.type, global_coh_elem.ghost_type)(ge);
    auto model_mat_index = coh_model.getMaterialByElement(
        global_coh_elem.type, global_coh_elem.ghost_type)(ge);
    const auto & material = coh_model.getMaterial(model_mat_index);

    // /// access cohesive material
    // const auto & mat = dynamic_cast<const MaterialCohesive &>(material);
    // const auto & normal_open_norm = mat.getNormalOpeningNorm(typecoh, gt);
    const auto normal_open_norm_it =
        make_view(material.getArray<Real>("normal_opening_norm", typecoh),
                  nb_quad_coh_elem)
            .begin();
    if (nb_quad_elem == 1) {
      /// coh el quad points are averaged to assign single opening to flow el
      Real aperture_av = 0.;
      for (UInt ip : arange(nb_quad_coh_elem)) {
        auto normal_open_norm_ip = normal_open_norm_it[le](ip);
        aperture_av += normal_open_norm_ip;
      }
      aperture_av /= nb_quad_coh_elem;
      auto & aperture_ip = aperture_it[ge](0);
      if (aperture_av < this->default_aperture) {
        aperture_ip = this->default_aperture;
      } else {
        aperture_ip = aperture_av;
      }
    } else {
      /// each flow el quad point gets aperture from corresponding coh el qp
      for (UInt ip : arange(nb_quad_coh_elem)) {
        auto normal_open_norm_ip = normal_open_norm_it[le](ip);
        auto & aperture_ip = aperture_it[ge](ip);
        if (normal_open_norm_ip < this->default_aperture) {
          aperture_ip = this->default_aperture;
        } else {
          aperture_ip = normal_open_norm_ip;
        }
      }
    }
  }

  // save previous aperture
  if (first_time && use_aperture_speed) {
    prev_aperture_on_qpoints.copy(aperture_on_qpoints);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::updateFluidElementsFromCohesive(
    const SolidMechanicsModelCohesive & coh_model) {
  AKANTU_DEBUG_IN();

  // get element types
  auto & coh_mesh = coh_model.getMesh();
  const auto dim = coh_mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto type = *coh_mesh.elementTypes(dim, gt, _ek_regular).begin();
  const ElementType type_facets = Mesh::getFacetType(type);
  const ElementType typecoh = FEEngine::getCohesiveElementType(type_facets);
  const auto nb_coh_elem = coh_mesh.getNbElement(typecoh);
  const auto nb_quad_coh_elem = coh_model.getFEEngine("CohesiveFEEngine")
                                    .getNbIntegrationPoints(typecoh, gt);

  auto type_fluid =
      *mesh.elementTypes(spatial_dimension, gt, _ek_regular).begin();
  const auto nb_fluid_elem = mesh.getNbElement(type_fluid);

  // check if there is need to add more fluid elements
  if (nb_coh_elem == nb_fluid_elem) {
    return;
  }

  NewElementsEvent new_facets;
  NewNodesEvent new_nodes;
  NewElementsEvent new_fluid_facets;

  /// loop on each cohesive element and check its damage
  for (auto && coh_el : arange(nb_coh_elem)) {

    // check if this cohesive element is present in fluid mesh
    auto & existing_coh_elements =
        mesh.getData<Element>("cohesive_elements", type_facets, gt);
    Element coh_element = {typecoh, coh_el, gt};
    auto it = existing_coh_elements.find(coh_element);
    if (it != UInt(-1)) {
      continue;
    }

    // check damage at the cohesive element
    // Element element{type_fluid, element_nb, gt};
    auto le = coh_model.getMaterialLocalNumbering(typecoh, gt)(coh_el);
    auto model_mat_index = coh_model.getMaterialByElement(typecoh, gt)(coh_el);
    const auto & material = coh_model.getMaterial(model_mat_index);
    const auto damage_it =
        make_view(material.getArray<Real>("damage", typecoh), nb_quad_coh_elem)
            .begin();
    Real damage_av = 0.;
    for (auto ip : arange(nb_quad_coh_elem)) {
      auto damage_ip = damage_it[le](ip);
      damage_av += damage_ip;
    }
    damage_av /= nb_quad_coh_elem;

    if (damage_av < this->insertion_damage) {
      continue;
    }

    // add fluid element to connectivity and corresponding nodes
    Array<UInt> & nodes_added = new_nodes.getList();
    auto & crack_facets = coh_mesh.getElementGroup("crack_facets");
    auto & nodes = mesh.getNodes();
    auto && coh_nodes = coh_mesh.getNodes();
    auto && coh_conn = coh_mesh.getConnectivity(typecoh, gt);
    Vector<UInt> conn(coh_conn.getNbComponent() / 2);
    auto coh_facet_conn_it = make_view(coh_mesh.getConnectivity(typecoh, gt),
                                       coh_conn.getNbComponent() / 2)
                                 .begin();
    auto cohesive_nodes_first_half = coh_facet_conn_it[2 * coh_el];
    for (auto c : akantu::arange(conn.size())) {
      auto coh_node = cohesive_nodes_first_half(c);
      Vector<Real> pos(mesh.getSpatialDimension());
      for (auto && data : enumerate(pos)) {
        std::get<1>(data) = coh_nodes(coh_node, std::get<0>(data));
      }
      auto idx = nodes.find(pos);
      if (idx == UInt(-1)) {
        nodes.push_back(pos);
        conn(c) = nodes.size() - 1;
        nodes_added.push_back(nodes.size() - 1);
      } else {
        conn(c) = idx;
      }
    }
    mesh.getData<Element>("cohesive_elements", type_facets, gt)
        .push_back(coh_element);
    mesh.getConnectivity(type_facets, gt).push_back(conn);
    Element fluid_facet{type_facets,
                        mesh.getConnectivity(type_facets, gt).size() - 1, gt};
    new_fluid_facets.getList().push_back(fluid_facet);

    // instantiate connectivity type for this facet type if not yet done
    if (not coh_mesh.getConnectivities().exists(type_facets, gt)) {
      coh_mesh.addConnectivityType(type_facets, gt);
    }

    // add two facets per cohesive element into cohesive mesh
    auto & facet_conn_mesh = coh_mesh.getConnectivity(type_facets, gt);
    for (auto f : akantu::arange(2)) {
      Vector<UInt> coh_nodes_one_side = coh_facet_conn_it[2 * coh_el + f];
      facet_conn_mesh.push_back(coh_nodes_one_side);
      Element new_facet{type_facets, facet_conn_mesh.size() - 1, gt};
      crack_facets.add(new_facet);
      new_facets.getList().push_back(new_facet);
    }
  }
  mesh.sendEvent(new_nodes);
  mesh.sendEvent(new_fluid_facets);
  MeshUtils::fillElementToSubElementsData(coh_mesh);
  coh_mesh.sendEvent(new_facets);

  AKANTU_DEBUG_OUT();
}

#endif

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::applyExternalFluxAtElementGroup(
    const Real & rate, const ElementGroup & source_facets,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && type :
       source_facets.elementTypes(spatial_dimension, ghost_type)) {
    const auto & element_ids = source_facets.getElements(type, ghost_type);

    const auto nb_fluid_elem = mesh.getNbElement(type);
    const auto type_size = element_ids.size();
    AKANTU_DEBUG_ASSERT(type_size <= nb_fluid_elem,
                        "Number of provided source facets exceeds total number "
                        "of flow elements");
    const auto fluid_conn_it =
        make_view(mesh.getConnectivity(type, ghost_type),
                  mesh.getConnectivity(type, ghost_type).getNbComponent())
            .begin();
    auto & source = getExternalFlux();
    source.clear();

    /// loop on each element of the group
    for (auto el : element_ids) {
      AKANTU_DEBUG_ASSERT(el <= nb_fluid_elem,
                          "Element number in the source facets group exceeds "
                          "maximum number of flow elements");

      const Vector<UInt> fluid_conn = fluid_conn_it[el];
      for (const auto & node : fluid_conn) {
        source(node) += rate / (fluid_conn.size() * element_ids.size());
      }
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Array<Real> FluidDiffusionModel::getQpointsCoord() {
  AKANTU_DEBUG_IN();

  const auto & nodes_coords = mesh.getNodes();
  const GhostType gt = _not_ghost;
  auto type_fluid =
      *mesh.elementTypes(spatial_dimension, gt, _ek_regular).begin();
  const auto nb_fluid_elem = mesh.getNbElement(type_fluid);
  auto & fem = this->getFEEngine();
  const auto nb_quad_points = fem.getNbIntegrationPoints(type_fluid, gt);

  Array<Real> quad_coords(nb_fluid_elem * nb_quad_points,
                          spatial_dimension + 1);
  fem.interpolateOnIntegrationPoints(nodes_coords, quad_coords,
                                     spatial_dimension + 1, type_fluid, gt);
  AKANTU_DEBUG_OUT();
  return quad_coords;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

std::shared_ptr<dumpers::Field> FluidDiffusionModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  auto field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> FluidDiffusionModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION(
        "Capacity lumped is a nodal field now stored in the DOF manager."
        "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["pressure"] = pressure.get();
  real_nodal_fields["pressure_rate"] = pressure_rate.get();
  real_nodal_fields["external_flux"] = external_flux.get();
  real_nodal_fields["internal_flux"] = internal_flux.get();
  // real_nodal_fields["increment"] = increment.get();

  std::shared_ptr<dumpers::Field> field =
      mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> FluidDiffusionModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*padding_flag*/, UInt /*spatial_dimension*/,
    ElementKind element_kind) {

  std::shared_ptr<dumpers::Field> field;

  if (field_name == "partitions") {
    field = mesh.createElementalField<UInt, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  } else if (field_name == "pressure_gradient") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(pressure_gradient);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        pressure_gradient, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  } else if (field_name == "permeability") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(permeability_on_qpoints);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        permeability_on_qpoints, group_name, this->spatial_dimension,
        element_kind, nb_data_per_elem);
  } else if (field_name == "aperture") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(aperture_on_qpoints);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        aperture_on_qpoints, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  } else if (field_name == "pressure_on_qpoints") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(pressure_on_qpoints);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        pressure_on_qpoints, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> FluidDiffusionModel::createElementalField(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const ElementKind & element_kind) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> FluidDiffusionModel::createNodalFieldBool(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> FluidDiffusionModel::createNodalFieldReal(
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

/* -------------------------------------------------------------------------
 */
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
  case SynchronizationTag::_fdm_pressure: {
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
inline void
FluidDiffusionModel::packData(CommunicationBuffer & buffer,
                              const Array<UInt> & indexes,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_fdm_pressure: {
      buffer << (*pressure)(index);
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
inline void FluidDiffusionModel::unpackData(CommunicationBuffer & buffer,
                                            const Array<UInt> & indexes,
                                            const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_fdm_pressure: {
      buffer >> (*pressure)(index);
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
  case SynchronizationTag::_fdm_pressure: {
    size += nb_nodes_per_element * sizeof(Real); // pressure
    break;
  }
  case SynchronizationTag::_fdm_gradient_pressure: {
    // pressure gradient
    size += getNbIntegrationPoints(elements) * sizeof(Real);
    size += nb_nodes_per_element * sizeof(Real); // nodal pressures
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
inline void
FluidDiffusionModel::packData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_fdm_pressure: {
    packNodalDataHelper(*pressure, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_fdm_gradient_pressure: {
    packElementalDataHelper(pressure_gradient, buffer, elements, true,
                            getFEEngine());
    packNodalDataHelper(*pressure, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void FluidDiffusionModel::unpackData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_fdm_pressure: {
    unpackNodalDataHelper(*pressure, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_fdm_gradient_pressure: {
    unpackElementalDataHelper(pressure_gradient, buffer, elements, true,
                              getFEEngine());
    unpackNodalDataHelper(*pressure, buffer, elements, mesh);

    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
void FluidDiffusionModel::injectIntoFacetsByCoord(const Vector<Real> & position,
                                                  const Real & injection_rate) {
  AKANTU_DEBUG_IN();

  auto dim = mesh.getSpatialDimension();
  const auto gt = _not_ghost;
  const auto & pos = mesh.getNodes();
  const auto pos_it = make_view(pos, dim).begin();
  auto & source = *external_flux;

  Vector<Real> bary_facet(dim);

  Real min_dist = std::numeric_limits<Real>::max();
  Element inj_facet;
  for_each_element(
      mesh,
      [&](auto && facet) {
        mesh.getBarycenter(facet, bary_facet);
        auto dist = std::abs(bary_facet.distance(position));
        if (dist < min_dist) {
          min_dist = dist;
          inj_facet = facet;
        }
      },
      _spatial_dimension = dim - 1);

  auto & facet_conn = mesh.getConnectivity(inj_facet.type, gt);
  auto nb_nodes_per_elem = facet_conn.getNbComponent();
  for (auto node : arange(nb_nodes_per_elem)) {
    source(facet_conn(inj_facet.element, node)) =
        injection_rate / nb_nodes_per_elem;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
