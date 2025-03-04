/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "structural_mechanics_model.hh"
#include "dof_manager.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "shape_structural.hh"
#include "sparse_matrix.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */
#include "dumpable_inline_impl.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
#include "group_manager_inline_impl.hh"
/* -------------------------------------------------------------------------- */
#include "structural_element_bernoulli_beam_2.hh"
#include "structural_element_bernoulli_beam_3.hh"
#include "structural_element_kirchhoff_shell.hh"
/* -------------------------------------------------------------------------- */
// #include "structural_mechanics_model_inline_impl.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt StructuralMechanicsModel::getNbDegreeOfFreedom(ElementType type) {
  return tuple_dispatch<ElementTypes_t<_ek_structural>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        return ElementClass<type>::getNbDegreeOfFreedom();
      },
      type);
}

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::StructuralMechanicsModel(Mesh & mesh, Int dim,
                                                   const ID & id)
    : Model(mesh, ModelType::_structural_mechanics_model, dim, id), f_m2a(1.0),
      stress("stress", id), element_material("element_material", id),
      set_ID("beam sets", id) {
  AKANTU_DEBUG_IN();

  registerFEEngineObject<MyFEEngineType>("StructuralMechanicsFEEngine", mesh,
                                         spatial_dimension);

  if (spatial_dimension == 2) {
    nb_degree_of_freedom = 3;
  } else if (spatial_dimension == 3) {
    nb_degree_of_freedom = 6;
  } else {
    AKANTU_TO_IMPLEMENT();
  }

  this->mesh.registerDumper<DumperParaview>("structural_mechanics_model", id,
                                            true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_structural);

  this->initDOFManager();

  this->dumper_default_element_kind = _ek_structural;

  mesh.getElementalData<Real>("extra_normal")
      .initialize(mesh, _element_kind = _ek_structural,
                  _nb_component = spatial_dimension, _with_nb_element = true,
                  _default_value = 0.);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::~StructuralMechanicsModel() = default;

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  // Initializing stresses
  ElementTypeMap<UInt> stress_components;

  for (auto && type : mesh.elementTypes(_spatial_dimension = _all_dimensions,
                      _element_kind = _ek_structural)) {
    // Getting number of components for each element type
    auto nb_components = tuple_dispatch<ElementTypes_t<_ek_structural>>(
        [&](auto && enum_type) {
          constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
          return ElementClass<type>::getNbStressComponents();
        },
        type);

    stress_components(nb_components, type);
  }

  stress.initialize(
      getFEEngine(), _spatial_dimension = _all_dimensions,
      _element_kind = _ek_structural,
      _nb_component = [&stress_components](ElementType type,
                                           GhostType /*unused*/) -> UInt {
        return stress_components(type);
      });
}

/* -------------------------------------------------------------------------- */

void StructuralMechanicsModel::initFEEngineBoundary() {
  /// this function should not be reimplemented
  /// we're just avoiding a call to Model::initFEEngineBoundary()
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::setTimeStep(Real time_step,
                                           const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);
  this->mesh.getDumper().setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initSolver(
    TimeStepSolverType time_step_solver_type, NonLinearSolverType /*unused*/) {
  AKANTU_DEBUG_IN();

  this->allocNodalField(displacement_rotation, nb_degree_of_freedom,
                        "displacement");
  this->allocNodalField(external_force, nb_degree_of_freedom, "external_force");
  this->allocNodalField(internal_force, nb_degree_of_freedom, "internal_force");
  this->allocNodalField(blocked_dofs, nb_degree_of_freedom, "blocked_dofs");

  auto & dof_manager = this->getDOFManager();

  if (not dof_manager.hasDOFs("displacement")) {
    dof_manager.registerDOFs("displacement", *displacement_rotation,
                             _dst_nodal);
    dof_manager.registerBlockedDOFs("displacement", *this->blocked_dofs);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(velocity, nb_degree_of_freedom, "velocity");
    this->allocNodalField(acceleration, nb_degree_of_freedom, "acceleration");

    if (!dof_manager.hasDOFsDerivatives("displacement", 1)) {
      dof_manager.registerDOFsDerivative("displacement", 1, *this->velocity);
      dof_manager.registerDOFsDerivative("displacement", 2,
                                         *this->acceleration);
    }

    /* Only allocate the mass if the "lumped" mode is ennabled.
     *  Also it is not a 1D array, but has an element for every
     *  DOF, which are most of the time equal, but makes handling
     *  some operations a bit simpler. */
    if (time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
      this->allocateLumpedMassArray();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initModel() {
  element_material.initialize(mesh, _element_kind = _ek_structural,
                              _default_value = 0, _with_nb_element = true);

  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  if (not need_to_reassemble_stiffness) {
    return;
  }

  if (not getDOFManager().hasMatrix("K")) {
    getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }

  this->getDOFManager().zeroMatrix("K");

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_structural)) {

    tuple_dispatch<ElementTypes_t<_ek_structural>>(
        [&](auto && enum_type) {
          constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
          this->assembleStiffnessMatrix<type>();
        },
        type);
  }

  need_to_reassemble_stiffness = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeStresses() {
  AKANTU_DEBUG_IN();

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_structural)) {
    tuple_dispatch<ElementTypes_t<_ek_structural>>(
        [&](auto && enum_type) {
          constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
          this->computeStressOnQuad<type>();
        },
        type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool StructuralMechanicsModel::allocateLumpedMassArray() {
  this->allocNodalField(this->mass, this->nb_degree_of_freedom, "lumped_mass");

  return true;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> StructuralMechanicsModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  return mesh.createNodalField(uint_nodal_fields[field_name], group_name);
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
StructuralMechanicsModel::createNodalFieldReal(const std::string & field_name,
                                               const std::string & group_name,
                                               bool padding_flag) {

  UInt n;
  if (spatial_dimension == 2) {
    n = 2;
  } else {
    n = 3;
  }

  UInt padding_size = 0;
  if (padding_flag) {
    padding_size = 3;
  }

  if (field_name == "displacement") {
    return mesh.createStridedNodalField(displacement_rotation.get(), group_name,
                                        n, 0, padding_size);
  }

  if (field_name == "velocity") {
    return mesh.createStridedNodalField(velocity.get(), group_name, n, 0,
                                        padding_size);
  }

  if (field_name == "acceleration") {
    return mesh.createStridedNodalField(acceleration.get(), group_name, n, 0,
                                        padding_size);
  }

  if (field_name == "rotation") {
    return mesh.createStridedNodalField(displacement_rotation.get(), group_name,
                                        nb_degree_of_freedom - n, n,
                                        padding_size);
  }

  if (field_name == "force") {
    return mesh.createStridedNodalField(external_force.get(), group_name, n, 0,
                                        padding_size);
  }
  if (field_name == "external_force") {
    return mesh.createStridedNodalField(external_force.get(), group_name, n, 0,
                                        padding_size);
  }

  if (field_name == "momentum") {
    return mesh.createStridedNodalField(external_force.get(), group_name,
                                        nb_degree_of_freedom - n, n,
                                        padding_size);
  }

  if (field_name == "internal_force") {
    return mesh.createStridedNodalField(internal_force.get(), group_name, n, 0,
                                        padding_size);
  }

  if (field_name == "internal_momentum") {
    return mesh.createStridedNodalField(internal_force.get(), group_name,
                                        nb_degree_of_freedom - n, n,
                                        padding_size);
  }

  if (field_name == "mass") {
    AKANTU_DEBUG_ASSERT(this->mass.get() != nullptr,
                        "The lumped mass matrix was not allocated.");
    return mesh.createStridedNodalField(this->mass.get(), group_name, n, 0,
                                        padding_size);
  }

  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> StructuralMechanicsModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*unused*/, Int spatial_dimension, ElementKind kind) {

  std::shared_ptr<dumpers::Field> field;

  if (field_name == "element_index_by_material") {
    field =
        mesh.createElementalField<Idx, Vector<Idx>, dumpers::ElementalField>(
            field_name, group_name, spatial_dimension, kind);
  }
  if (field_name == "stress") {
    ElementTypeMap<Int> nb_data_per_elem = this->mesh.getNbDataPerElem(stress);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        stress, group_name, this->spatial_dimension, kind, nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
/* Virtual methods from SolverCallback */
/* -------------------------------------------------------------------------- */
/// get the type of matrix needed
MatrixType StructuralMechanicsModel::getMatrixType(const ID & /*id*/) const {
  return _symmetric;
}

/// callback to assemble a Matrix
void StructuralMechanicsModel::assembleMatrix(const ID & id) {
  if (id == "K") {
    assembleStiffnessMatrix();
  } else if (id == "M") {
    assembleMassMatrix();
  }
}

/// callback to assemble a lumped Matrix
void StructuralMechanicsModel::assembleLumpedMatrix(const ID & id) {
  if ("M" == id) {
    this->assembleLumpedMassMatrix();
  }
  return;
}

/// callback to assemble the residual StructuralMechanicsModel::(rhs)
void StructuralMechanicsModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  auto & dof_manager = getDOFManager();

  assembleInternalForce();

  /* This is essentially a summing up of forces
   * first the external forces are counted for and then stored inside the
   * residual.
   */
  dof_manager.assembleToResidual("displacement", *external_force, 1);
  dof_manager.assembleToResidual("displacement", *internal_force, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_IN();

  auto & dof_manager = this->getDOFManager();

  if ("external" == residual_part) {
    dof_manager.assembleToResidual("displacement", *this->external_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  if ("internal" == residual_part) {
    this->assembleInternalForce();
    dof_manager.assembleToResidual("displacement", *this->internal_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_CUSTOM_EXCEPTION(
      debug::SolverCallbackResidualPartUnknown(residual_part));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Virtual methods from Model                                                 */
/* -------------------------------------------------------------------------- */
/// get some default values for derived classes
std::tuple<ID, TimeStepSolverType>
StructuralMechanicsModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  case _implicit_dynamic: {
    return std::make_tuple("implicit", TimeStepSolverType::_dynamic);
  }
  case _explicit_lumped_mass: { // Taken from the solid mechanics part
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _explicit_consistent_mass: { // Taken from the solid mechanics part
    return std::make_tuple("explicit", TimeStepSolverType::_dynamic);
  }
  default:
    std::cout << "UNKOWN." << std::endl;
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }
}

/* ------------------------------------------------------------------------ */
ModelSolverOptions StructuralMechanicsModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: { // Taken from the solid mechanic
                                              // part
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type =
        NonLinearSolverType::_newton_raphson; // _linear;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
#if 1
  case TimeStepSolverType::_dynamic: { // Copied from solid
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["displacement"] =
          IntegrationSchemeType::_central_difference;
      options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["displacement"] =
          IntegrationSchemeType::_trapezoidal_rule_2;
      options.solution_type["displacement"] = IntegrationScheme::_displacement;
    }
    break;
  }
#else
  case TimeStepSolverType::_dynamic: { // Original
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_trapezoidal_rule_2;
    options.solution_type["displacement"] = IntegrationScheme::_displacement;
    break;
  }
#endif
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleInternalForce() {

  internal_force->zero();
  computeStresses();

  for (auto type : mesh.elementTypes(_spatial_dimension = _all_dimensions,
                   _element_kind = _ek_structural)) {
    assembleInternalForce(type, _not_ghost);
    // assembleInternalForce(type, _ghost);
  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleInternalForce(ElementType type,
                                                     GhostType gt) {
  auto & fem = getFEEngine();
  auto & sigma = stress(type, gt);
  auto ndof = getNbDegreeOfFreedom(type);
  auto nb_nodes = mesh.getNbNodesPerElement(type);
  auto ndof_per_elem = ndof * nb_nodes;

  Array<Real> BtSigma(fem.getNbIntegrationPoints(type) *
                          mesh.getNbElement(type),
                      ndof_per_elem, "BtSigma");
  fem.computeBtD(sigma, BtSigma, type, gt);

  Array<Real> intBtSigma(0, ndof_per_elem, "intBtSigma");
  fem.integrate(BtSigma, intBtSigma, ndof_per_elem, type, gt);

  getDOFManager().assembleElementalArrayLocalArray(intBtSigma, *internal_force,
                                                   type, gt, -1.);
}

/* -------------------------------------------------------------------------- */
Real StructuralMechanicsModel::getKineticEnergy() {

  const UInt nb_nodes = mesh.getNbNodes();
  const UInt nb_degree_of_freedom = this->nb_degree_of_freedom;
  Real ekin = 0.; // used to sum up energy (is divided by two at the very end)

  if (this->getDOFManager().hasLumpedMatrix("M")) {
    /* This code computes the kinetic energy for the case when the mass is
     * lumped. It is based on the solid mechanic equivalent.
     */
    AKANTU_DEBUG_ASSERT(this->mass != nullptr,
                        "The lumped mass is not allocated.");

    this->assembleLumpedMatrix("M");

    /* Iterating over all nodes.
     *   Important the velocity and mass also contains the rotational parts.
     *   However, they can be handled in an uniform way. */
    for (auto && data :
         zip(arange(nb_nodes), make_view(*this->velocity, nb_degree_of_freedom),
             make_view(*this->mass, nb_degree_of_freedom))) {
      const UInt n = std::get<0>(data); // This is the ID of the current node

      // Only handle the node if it belongs to us.
      if (not mesh.isLocalOrMasterNode(n)) {
        continue;
      }

      const auto & velocity = std::get<1>(data);
      const auto & mass = std::get<2>(data);
      Real mv2 = 0.;

      for (auto && node_data : zip(velocity, mass)) {
        mv2 += Math::pow<2>(std::get<0>(node_data)) * std::get<1>(node_data);
      }

      ekin += mv2;
    }
  } else if (this->getDOFManager().hasMatrix("M")) {
    /* Handle the case where no lumped mass is there.
     * This is basically the original code.
     */
    this->assembleMassMatrix();

    Array<Real> Mv(nb_nodes, nb_degree_of_freedom);
    this->getDOFManager().assembleMatMulVectToArray("displacement", "M",
                                                    *this->velocity, Mv);

    for (auto && data :
         zip(arange(nb_nodes), make_view(Mv, nb_degree_of_freedom),
             make_view(*this->velocity, nb_degree_of_freedom))) {
      // only consider the node if we are belonging to it
      if (mesh.isLocalOrMasterNode(std::get<0>(data))) {
        ekin += std::get<2>(data).dot(std::get<1>(data));
      }
    }
  } else {
    /* This is the case where no mass is present, for whatever reason, such as
     * the static case. We handle it specially be returning directly zero.
     * However, by doing that there will not be a synchronizing event as in the
     * other cases. Which is faster, but could be a problem in case the user
     * expects this.
     *
     * Another not is, that the solid mechanics part, would generate an error in
     * this clause. But, since the original implementation of the structural
     * part, did not do that, I, Philip, decided to refrain from that. However,
     * it is an option that should be considered.
     */
    return 0.;
  }

  // Sum up across the comunicator
  mesh.getCommunicator().allReduce(ekin, SynchronizerOperation::_sum);

  return ekin / 2.; // finally divide the energy by two
}

/* -------------------------------------------------------------------------- */
Real StructuralMechanicsModel::getPotentialEnergy() {
  Real epot = 0.;
  UInt nb_nodes = mesh.getNbNodes();

  // if stiffness matrix is not assembled, do it
  this->assembleStiffnessMatrix();

  Array<Real> Ku(nb_nodes, nb_degree_of_freedom);
  this->getDOFManager().assembleMatMulVectToArray(
      "displacement", "K", *this->displacement_rotation, Ku);

  for (auto && data :
       zip(arange(nb_nodes), make_view(Ku, nb_degree_of_freedom),
           make_view(*this->displacement_rotation, nb_degree_of_freedom))) {
    epot += std::get<2>(data).dot(std::get<1>(data)) *
            static_cast<Real>(mesh.isLocalOrMasterNode(std::get<0>(data)));
  }

  mesh.getCommunicator().allReduce(epot, SynchronizerOperation::_sum);

  return epot / 2.;
}

/* -------------------------------------------------------------------------- */
Real StructuralMechanicsModel::getEnergy(const ID & energy) {
  if (energy == "kinetic") {
    return getKineticEnergy();
  }

  if (energy == "potential") {
    return getPotentialEnergy();
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeForcesByLocalTractionArray(
    const Array<Real> & tractions, ElementType type) {
  AKANTU_DEBUG_IN();

  auto nb_element = getFEEngine().getMesh().getNbElement(type);
  auto nb_nodes_per_element =
      getFEEngine().getMesh().getNbNodesPerElement(type);
  auto nb_quad = getFEEngine().getNbIntegrationPoints(type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(
      Mesh::getSpatialDimension(type) == getFEEngine().getElementDimension(),
      "element type dimension does not match the dimension of boundaries : "
          << getFEEngine().getElementDimension()
          << " != " << Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(
      tractions.size() == nb_quad * nb_element,
      "the size of the vector should be the total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(tractions.getNbComponent() == nb_degree_of_freedom,
                      "the number of components should be the spatial "
                      "dimension of the problem");

  Array<Real> Ntbs(nb_element * nb_quad,
                   nb_degree_of_freedom * nb_nodes_per_element);

  auto & fem = getFEEngine();
  fem.computeNtb(tractions, Ntbs, type);

  // allocate the vector that will contain the integrated values
  auto name = id + std::to_string(type) + ":integral_boundary";
  Array<Real> int_funct(nb_element, nb_degree_of_freedom * nb_nodes_per_element,
                        name);

  // do the integration
  getFEEngine().integrate(Ntbs, int_funct,
                          nb_degree_of_freedom * nb_nodes_per_element, type);

  // assemble the result into force vector
  getDOFManager().assembleElementalArrayLocalArray(int_funct, *external_force,
                                                   type, _not_ghost, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeForcesByGlobalTractionArray(
    const Array<Real> & traction_global, ElementType type) {
  AKANTU_DEBUG_IN();

  auto nb_element = mesh.getNbElement(type);
  auto nb_quad = getFEEngine().getNbIntegrationPoints(type);

  Array<Real> traction_local(nb_element * nb_quad, nb_degree_of_freedom,
                             id + ":structuralmechanics:imposed_linear_load");

  auto R_it = getFEEngineClass<MyFEEngineType>()
                  .getShapeFunctions()
                  .getRotations(type)
                  .begin(nb_degree_of_freedom, nb_degree_of_freedom);

  auto Te_it = traction_global.begin(nb_degree_of_freedom);
  auto te_it = traction_local.begin(nb_degree_of_freedom);

  for (Int e = 0; e < nb_element; ++e, ++R_it) {
    for (Int q = 0; q < nb_quad; ++q, ++Te_it, ++te_it) {
      // turn the traction in the local referential
      *te_it = *R_it * *Te_it;
    }
  }

  computeForcesByLocalTractionArray(traction_local, type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::afterSolveStep(bool converged) {
  if (converged) {
    assembleInternalForce();
  }
}

} // namespace akantu
