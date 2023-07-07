/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "phase_field_model.hh"
#include "aka_common.hh"
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
#include <utility>
/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseFieldModel::PhaseFieldModel(Mesh & mesh, Int dim, const ID & id,
                                 std::shared_ptr<DOFManager> dof_manager,
                                 ModelType model_type)
    : Model(mesh, model_type, dim, id),
      phasefield_index("phasefield index", id),
      phasefield_local_numbering("phasefield local numbering", id) {

  AKANTU_DEBUG_IN();

  this->initDOFManager(std::move(dof_manager));

  this->registerFEEngineObject<FEEngineType>("PhaseFieldFEEngine", mesh,
                                             Model::spatial_dimension);

  this->mesh.registerDumper<DumperParaview>("phase_field", id, true);
  this->mesh.addDumpMesh(mesh, Model::spatial_dimension, _not_ghost,
                         _ek_regular);

  phasefield_selector =
      std::make_shared<DefaultPhaseFieldSelector>(phasefield_index);

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_phasefield_id);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_pfm_damage);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_for_dump);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PhaseFieldModel::~PhaseFieldModel() = default;

/* -------------------------------------------------------------------------- */
MatrixType PhaseFieldModel::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K") {
    return _symmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initFullImpl(const ModelOptions & options) {
  phasefield_index.initialize(mesh, _element_kind = _ek_not_defined,
                              _default_value = Idx(-1),
                              _with_nb_element = true);
  phasefield_local_numbering.initialize(mesh, _element_kind = _ek_not_defined,
                                        _with_nb_element = true);

  Model::initFullImpl(options);

  // initialize the phasefields
  if (!this->parser.getLastParsedFile().empty()) {
    this->instantiatePhaseFields();
    this->initPhaseFields();
  }

  this->initBC(*this, *damage, *external_force);
}

/* -------------------------------------------------------------------------- */
PhaseField &
PhaseFieldModel::registerNewPhaseField(const ParserSection & section) {
  std::string phase_name;
  std::string phase_type = section.getName();
  std::string opt_param = section.getOption();

  try {
    std::string tmp = section.getParameter("name");
    phase_name = tmp; /** this can seam weird, but there is an ambiguous
                       * operator overload that i couldn't solve. @todo remove
                       * the weirdness of this code
                       */
  } catch (debug::Exception &) {
    AKANTU_ERROR("A phasefield of type \'"
                 << phase_type
                 << "\' in the input file has been defined without a name!");
  }
  PhaseField & phase =
      this->registerNewPhaseField(phase_name, phase_type, opt_param);

  phase.parseSection(section);

  return phase;
}

/* -------------------------------------------------------------------------- */
PhaseField & PhaseFieldModel::registerNewPhaseField(const ID & phase_name,
                                                    const ID & phase_type,
                                                    const ID & opt_param) {
  AKANTU_DEBUG_ASSERT(phasefields_names_to_id.find(phase_name) ==
                          phasefields_names_to_id.end(),
                      "A phasefield with this name '"
                          << phase_name << "' has already been registered. "
                          << "Please use unique names for phasefields");

  Int phase_count = phasefields.size();
  phasefields_names_to_id[phase_name] = phase_count;

  std::stringstream sstr_phase;
  sstr_phase << this->id << ":" << phase_count << ":" << phase_type;
  ID mat_id = sstr_phase.str();

  std::unique_ptr<PhaseField> phase = PhaseFieldFactory::getInstance().allocate(
      phase_type, opt_param, *this, mat_id);

  phasefields.push_back(std::move(phase));

  return *(phasefields.back());
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::instantiatePhaseFields() {

  ParserSection model_section;
  bool is_empty;
  std::tie(model_section, is_empty) = this->getParserSection();

  if (not is_empty) {
    auto model_phasefields =
        model_section.getSubSections(ParserType::_phasefield);
    for (const auto & section : model_phasefields) {
      this->registerNewPhaseField(section);
    }
  }

  auto sub_sections = this->parser.getSubSections(ParserType::_phasefield);
  for (const auto & section : sub_sections) {
    this->registerNewPhaseField(section);
  }

  if (phasefields.empty()) {
    AKANTU_EXCEPTION("No phasefields where instantiated for the model"
                     << getID());
  }
  are_phasefields_instantiated = true;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initPhaseFields() {
  AKANTU_DEBUG_ASSERT(phasefields.size() != 0, "No phasefield to initialize !");

  if (!are_phasefields_instantiated) {
    instantiatePhaseFields();
  }

  this->assignPhaseFieldToElements();

  for (auto & phasefield : phasefields) {
    /// init internals properties
    phasefield->initPhaseField();
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assignPhaseFieldToElements(
    const ElementTypeMapArray<Idx> * filter) {

  for_each_element(
      mesh,
      [&](auto && element) {
        Int phase_index = (*phasefield_selector)(element);
        AKANTU_DEBUG_ASSERT(
            phase_index < Int(phasefields.size()),
            "The phasefield selector returned an index that does not exists");
        phasefield_index(element) = phase_index;
      },
      _element_filter = filter, _ghost_type = _not_ghost);

  for_each_element(
      mesh,
      [&](auto && element) {
        auto phase_index = phasefield_index(element);
        auto index = phasefields[phase_index]->addElement(element);
        phasefield_local_numbering(element) = index;
      },
      _element_filter = filter, _ghost_type = _not_ghost);

  // synchronize the element phasefield arrays
  this->synchronize(SynchronizationTag::_phasefield_id);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleMatrix(const ID & matrix_id) {

  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else {
    AKANTU_ERROR("Unknown Matrix ID for PhaseFieldModel : " << matrix_id);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::predictor() {
  // AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::corrector() {
  // AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initSolver(TimeStepSolverType time_step_solver_type,
                                 NonLinearSolverID /*unused*/) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->damage, 1, "damage");
  this->allocNodalField(this->external_force, 1, "external_force");
  this->allocNodalField(this->internal_force, 1, "internal_force");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");
  this->allocNodalField(this->previous_damage, 1, "previous_damage");

  if (!dof_manager.hasDOFs("damage")) {
    dof_manager.registerDOFs("damage", *this->damage, _dst_nodal);
    dof_manager.registerBlockedDOFs("damage", *this->blocked_dofs);
    dof_manager.registerDOFsPrevious("damage", *this->previous_damage);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic) {
    AKANTU_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
FEEngine & PhaseFieldModel::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(getFEEngineClassBoundary<FEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType PhaseFieldModel::getDefaultSolverType() const {
  return TimeStepSolverType::_static;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
PhaseFieldModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _explicit_consistent_mass: {
    return std::make_tuple("explicit", TimeStepSolverType::_dynamic);
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
ModelSolverOptions PhaseFieldModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = "lumped";
    options.integration_scheme_type["damage"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["damage"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = "newton_raphson";
    options.integration_scheme_type["damage"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["damage"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    options.non_linear_solver_type = "newton_raphson";
    options.integration_scheme_type["damage"] =
        IntegrationSchemeType::_backward_euler;
    options.solution_type["damage"] = IntegrationScheme::_damage;
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
Real PhaseFieldModel::getEnergy() {
  AKANTU_DEBUG_IN();

  Real energy = 0.;
  for (auto & phasefield : phasefields) {
    energy += phasefield->getEnergy();
  }

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real PhaseFieldModel::getEnergy(ElementType type, Idx index) {
  AKANTU_DEBUG_IN();

  Idx phase_index = this->phasefield_index(type, _not_ghost)(index);
  Idx phase_loc_num = this->phasefield_local_numbering(type, _not_ghost)(index);
  Real energy = this->phasefields[phase_index]->getEnergy(
      Element{type, phase_loc_num, _not_ghost});

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real PhaseFieldModel::getEnergy(const ID & group_id) {
  auto && group = mesh.getElementGroup(group_id);
  auto energy = 0.;
  for (auto && type : group.elementTypes()) {
    for (auto el : group.getElementsIterable(type)) {
      energy += getEnergy(el);
    }
  }

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  return energy;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::beforeSolveStep() {
  for (auto & phasefield : phasefields) {
    phasefield->beforeSolveStep();
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::afterSolveStep(bool converged) {
  if (not converged) {
    return;
  }

  for (auto && values : zip(*damage, *previous_damage)) {
    auto & dam = std::get<0>(values);
    auto & prev_dam = std::get<1>(values);

    prev_dam = dam;
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }

  this->getDOFManager().zeroMatrix("K");

  for (auto & phasefield : phasefields) {
    phasefield->assembleStiffnessMatrix(_not_ghost);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleResidual() {
  this->assembleInternalForces();

  this->getDOFManager().assembleToResidual("damage", *this->external_force, 1);
  this->getDOFManager().assembleToResidual("damage", *this->internal_force, 1);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleInternalForces() {
  AKANTU_DEBUG_INFO("Assemble the internal forces");

  this->internal_force->zero();

  this->synchronize(SynchronizationTag::_pfm_damage);

  for (auto & phasefield : phasefields) {
    phasefield->computeAllDrivingForces(_not_ghost);
  }

  // assemble the forces due to local driving forces
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for (auto & phasefield : phasefields) {
    phasefield->assembleInternalForces(_not_ghost);
  }

  // assemble the forces due to local driving forces
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for (auto & phasefield : phasefields) {
    phasefield->assembleInternalForces(_ghost);
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleLumpedMatrix(const ID & /*matrix_id*/) {}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

  this->mesh.getDumper("phase_field").setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
Int PhaseFieldModel::getNbData(const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  Int size = 0;
  Int nb_nodes_per_element = 0;

  for (const Element & el : elements) {
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case SynchronizationTag::_phasefield_id: {
    size += elements.size() * sizeof(Int);
    break;
  }
  case SynchronizationTag::_for_dump: {
    // damage
    size += nb_nodes_per_element * sizeof(Real);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    size += nb_nodes_per_element * sizeof(Real);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  return size;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::packData(CommunicationBuffer & buffer,
                               const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_phasefield_id: {
    packElementalDataHelper(phasefield_index, buffer, elements, false,
                            getFEEngine());
    break;
  }
  case SynchronizationTag::_for_dump: {
    packNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    packNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::unpackData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_phasefield_id: {
    for (auto && element : elements) {
      Idx recv_phase_index;
      buffer >> recv_phase_index;
      Idx & phase_index = phasefield_index(element);
      if (phase_index != Idx(-1)) {
        continue;
      }

      // add ghosts element to the correct phasefield
      phase_index = recv_phase_index;
      Idx index = phasefields[phase_index]->addElement(element);
      phasefield_local_numbering(element) = index;
    }
    break;
  }
  case SynchronizationTag::_for_dump: {
    unpackNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    unpackNodalDataHelper(*damage, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Int PhaseFieldModel::getNbData(const Array<Idx> & indexes,
                               const SynchronizationTag & tag) const {
  Int size = 0;
  Int nb_nodes = indexes.size();

  switch (tag) {
  case SynchronizationTag::_for_dump: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
  return size;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::packData(CommunicationBuffer & buffer,
                               const Array<Idx> & indexes,
                               const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_for_dump: {
    packDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    packDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::unpackData(CommunicationBuffer & buffer,
                                 const Array<Idx> & indexes,
                                 const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_for_dump: {
    unpackDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  case SynchronizationTag::_pfm_damage: {
    unpackDOFDataHelper(*damage, buffer, indexes);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
PhaseFieldModel::createNodalFieldBool(const std::string & field_name,
                                      const std::string & group_name,
                                      bool /*unused*/) {
  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  return mesh.createNodalField(uint_nodal_fields[field_name], group_name);

  std::shared_ptr<dumpers::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
PhaseFieldModel::createNodalFieldReal(const std::string & field_name,
                                      const std::string & group_name,
                                      bool /*unused*/) {
  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["damage"] = damage.get();
  real_nodal_fields["external_force"] = external_force.get();
  real_nodal_fields["internal_force"] = internal_force.get();

  return mesh.createNodalField(real_nodal_fields[field_name], group_name);

  std::shared_ptr<dumpers::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PhaseFieldModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*unused*/, Int /*unused*/, ElementKind element_kind) {

  if (field_name == "partitions") {
    return mesh.createElementalField<Int, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  }

  std::shared_ptr<dumpers::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
ElementTypeMapArray<Real> &
PhaseFieldModel::flattenInternal(const std::string & field_name,
                                 ElementKind kind, const GhostType ghost_type) {
  auto key = std::make_pair(field_name, kind);

  ElementTypeMapArray<Real> * internal_flat;

  auto it = this->registered_internals.find(key);
  if (it == this->registered_internals.end()) {
    auto internal =
        std::make_unique<ElementTypeMapArray<Real>>(field_name, this->id);

    internal_flat = internal.get();
    this->registered_internals[key] = std::move(internal);
  } else {
    internal_flat = it->second.get();
  }

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, kind)) {
    if (internal_flat->exists(type, ghost_type)) {
      auto & internal = (*internal_flat)(type, ghost_type);
      internal.resize(0);
    }
  }

  for (auto & phasefield : phasefields) {
    if (phasefield->isInternal<Real>(field_name, kind)) {
      phasefield->flattenInternal(field_name, *internal_flat, ghost_type, kind);
    }
  }

  return *internal_flat;
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::inflateInternal(const std::string & field_name,
                                      const ElementTypeMapArray<Real> & field,
                                      ElementKind kind, GhostType ghost_type) {

  for (auto & phasefield : phasefields) {
    if (phasefield->isInternal<Real>(field_name, kind)) {
      phasefield->inflateInternal(field_name, field, ghost_type, kind);
    } else {
      AKANTU_ERROR("A internal of name \'"
                   << field_name
                   << "\' has not been defined in the phasefield");
    }
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "Phase Field Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << Model::spatial_dimension
         << std::endl;
  stream << space << " + fem [" << std::endl;
  getFEEngine().printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;
  stream << space << " + nodals information [" << std::endl;
  damage->printself(stream, indent + 2);
  external_force->printself(stream, indent + 2);
  internal_force->printself(stream, indent + 2);
  blocked_dofs->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + phasefield information [" << std::endl;
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

} // namespace akantu
