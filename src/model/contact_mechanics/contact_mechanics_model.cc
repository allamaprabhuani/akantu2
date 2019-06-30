/**
 * @file   coontact_mechanics_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Contact mechanics model
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
#include "contact_mechanics_model.hh"
#include "boundary_condition_functor.hh"
#include "dumpable_inline_impl.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ContactMechanicsModel::ContactMechanicsModel(
    Mesh & mesh, UInt dim, const ID & id, const MemoryID & memory_id,
    std::shared_ptr<DOFManager> dof_manager, const ModelType model_type)
    : Model(mesh, model_type, dof_manager, dim, id, memory_id) {

  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<MyFEEngineType>("ContactMechanicsModel", mesh,
                                               Model::spatial_dimension);
#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("contact_mechanics", id, true);
  this->mesh.addDumpMeshToDumper("contact_mechanics", mesh,
                                 Model::spatial_dimension - 1, _not_ghost,
                                 _ek_regular);
#endif

  this->registerDataAccessor(*this);

  this->detector =
      std::make_unique<ContactDetector>(this->mesh, id + ":contact_detector");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactMechanicsModel::~ContactMechanicsModel() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::initFullImpl(const ModelOptions & options) {

  Model::initFullImpl(options);

  // initalize the resolutions
  if (this->parser.getLastParsedFile() != "") {
    this->instantiateResolutions();
  }

  this->initResolutions();

  this->initBC(*this, *displacement, *displacement_increment, *external_force);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::instantiateResolutions() {
  ParserSection model_section;
  bool is_empty;
  std::tie(model_section, is_empty) = this->getParserSection();

  if (not is_empty) {
    auto model_resolutions =
        model_section.getSubSections(ParserType::_contact_resolution);
    for (const auto & section : model_resolutions) {
      this->registerNewResolution(section);
    }
  }

  auto sub_sections =
      this->parser.getSubSections(ParserType::_contact_resolution);
  for (const auto & section : sub_sections) {
    this->registerNewResolution(section);
  }

  if (resolutions.empty())
    AKANTU_EXCEPTION("No contact resolutions where instantiated for the model"
                     << getID());
  are_resolutions_instantiated = true;
}

/* -------------------------------------------------------------------------- */
Resolution &
ContactMechanicsModel::registerNewResolution(const ParserSection & section) {
  std::string res_name;
  std::string res_type = section.getName();
  std::string opt_param = section.getOption();

  try {
    std::string tmp = section.getParameter("name");
    res_name = tmp; /** this can seem weird, but there is an ambiguous operator
                     * overload that i couldn't solve. @todo remove the
                     * weirdness of this code
                     */
  } catch (debug::Exception &) {
    AKANTU_ERROR("A contact resolution of type \'"
                 << res_type
                 << "\' in the input file has been defined without a name!");
  }
  Resolution & res = this->registerNewResolution(res_name, res_type, opt_param);

  res.parseSection(section);

  return res;
}

/* -------------------------------------------------------------------------- */
Resolution & ContactMechanicsModel::registerNewResolution(
    const ID & res_name, const ID & res_type, const ID & opt_param) {
  AKANTU_DEBUG_ASSERT(resolutions_names_to_id.find(res_name) ==
                          resolutions_names_to_id.end(),
                      "A resolution with this name '"
                          << res_name << "' has already been registered. "
                          << "Please use unique names for resolutions");

  UInt res_count = resolutions.size();
  resolutions_names_to_id[res_name] = res_count;

  std::stringstream sstr_res;
  sstr_res << this->id << ":" << res_count << ":" << res_type;
  ID res_id = sstr_res.str();

  std::unique_ptr<Resolution> resolution =
      ResolutionFactory::getInstance().allocate(res_type, spatial_dimension,
                                                opt_param, *this, res_id);

  resolutions.push_back(std::move(resolution));

  return *(resolutions.back());
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::initResolutions() {
  AKANTU_DEBUG_ASSERT(resolutions.size() != 0,
                      "No resolutions to initialize !");

  if (!are_resolutions_instantiated)
    instantiateResolutions();

  // \TODO check if each resolution needs a initResolution() method
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::initModel() {
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
FEEngine & ContactMechanicsModel::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(
      getFEEngineClassBoundary<MyFEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::initSolver(
    TimeStepSolverType /*time_step_solver_type*/, NonLinearSolverType) {

  // for alloc type of solvers
  this->allocNodalField(this->displacement, spatial_dimension, "displacement");
  this->allocNodalField(this->displacement_increment, spatial_dimension,
                        "displacement_increment");
  this->allocNodalField(this->internal_force, spatial_dimension,
                        "internal_force");
  this->allocNodalField(this->external_force, spatial_dimension,
                        "external_force");

  this->allocNodalField(this->gaps, 1, "gaps");
  this->allocNodalField(this->previous_gaps, 1, "previous_gaps");
  this->allocNodalField(this->nodal_area, 1, "areas");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");

  this->allocNodalField(this->normals, spatial_dimension, "normals");
  this->allocNodalField(this->tangents, spatial_dimension, "tangents");

  // todo register multipliers as dofs for lagrange multipliers
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
ContactMechanicsModel::getDefaultSolverID(const AnalysisMethod & method) {

  switch (method) {
  case _explicit_contact: {
    return std::make_tuple("explicit_contact", TimeStepSolverType::_static);
  }
  case _implicit_contact: {
    return std::make_tuple("implicit_contact", TimeStepSolverType::_static);
  }
  case _explicit_dynamic_contact: {
    return std::make_tuple("explicit_dynamic_contact",
                           TimeStepSolverType::_dynamic_lumped);
    break;
  }
  default:
    return std::make_tuple("unkown", TimeStepSolverType::_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions ContactMechanicsModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type =
        NonLinearSolverType::_newton_raphson_contact;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
    break;
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  /* ------------------------------------------------------------------------ */
  // computes the internal forces
  this->assembleInternalForces();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->internal_force, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the contact forces");

  UInt nb_nodes = mesh.getNbNodes();
  this->internal_force->clear();
  internal_force->resize(nb_nodes, 0.);

  // assemble the forces due to contact
  auto assemble = [&](auto && ghost_type) {
    for (auto & resolution : resolutions) {
      resolution->assembleInternalForces(ghost_type);
    }
  };

  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  assemble(_not_ghost);

  // assemble the stresses due to ghost elements
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  assemble(_ghost);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::search() {

  this->detector->search(this->contact_map);

  for (auto & entry : contact_map) {
    auto & element = entry.second;

    if (element.gap < 0)
      element.gap = std::abs(element.gap);
    else
      element.gap = -element.gap;
  }

  this->assembleFieldsFromContactMap();

  this->computeNodalAreas();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::search(Array<Real> & increment) {

  this->contact_map.clear();

  this->detector->search(this->contact_map);

  for (auto & entry : contact_map) {

    auto & element = entry.second;
    const auto & connectivity = element.connectivity;
    const auto & normal = element.normal;

    auto slave_node = connectivity[0];
    Vector<Real> u_slave(spatial_dimension);

    auto master_node = connectivity[1];
    Vector<Real> u_master(spatial_dimension);

    for (UInt s : arange(spatial_dimension)) {
      u_slave(s) = increment(slave_node, s);
      u_master(s) = increment(master_node, s);
    }

    // todo check this initial guess
    auto u = (u_master.norm() > u_slave.norm()) ? u_master : u_slave;
    Real uv = Math::vectorDot(u.storage(), normal.storage(), spatial_dimension);

    if (element.gap - uv <= 0) {
      element.gap = std::abs(element.gap - uv);
    } else {
      element.gap = 0.0;
    }
  }

  this->assembleFieldsFromContactMap();

  this->computeNodalAreas();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleFieldsFromContactMap() {

  UInt nb_nodes = mesh.getNbNodes();

  this->previous_gaps->copy(*(this->gaps));
  this->gaps->clear();

  gaps->resize(nb_nodes, 0.);
  normals->resize(nb_nodes, 0.);
  tangents->resize(nb_nodes, 0.);

  if (this->contact_map.empty())
    return;

  for (auto & entry : contact_map) {
    const auto & element = entry.second;
    auto connectivity = element.connectivity;
    auto node = connectivity(0);

    (*gaps)[node] = element.gap;

    for (UInt i = 0; i < spatial_dimension; ++i) {
      (*normals)(node, i) = element.normal[i];
      (*tangents)(node, i) = element.tangents(0, i);
    }
  }
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::computeNodalAreas() {

  UInt nb_nodes = mesh.getNbNodes();

  this->nodal_area->clear();
  this->external_force->clear();

  nodal_area->resize(nb_nodes, 0.);
  external_force->resize(nb_nodes, 0.);

  auto & fem_boundary = this->getFEEngineBoundary(); 
  fem_boundary.computeNormalsOnIntegrationPoints(_not_ghost);
  fem_boundary.computeNormalsOnIntegrationPoints(_ghost);
  
  this->applyBC(
      BC::Neumann::FromHigherDim(Matrix<Real>::eye(spatial_dimension, 1)),
      mesh.getElementGroup("contact_surface"));

  for (auto && tuple :
       zip(*nodal_area, make_view(*external_force, spatial_dimension))) {
    auto & area = std::get<0>(tuple);
    auto & force = std::get<1>(tuple);

    for (auto & f : force)
      area += pow(f, 2);

    area = sqrt(area);
  }

  this->external_force->clear();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::beforeSolveStep() {}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::afterSolveStep() {}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "Contact Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << Model::spatial_dimension
         << std::endl;
  stream << space << " + fem [" << std::endl;
  getFEEngine().printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + resolutions [" << std::endl;
  for (auto & resolution : resolutions) {
    resolution->printself(stream, indent + 1);
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
MatrixType ContactMechanicsModel::getMatrixType(const ID & matrix_id) {

  if (matrix_id == "K")
    return _symmetric;

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  }
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }

  for (auto & resolution : resolutions) {
    resolution->assembleStiffnessMatrix(_not_ghost);
  }
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleLumpedMatrix(const ID & /*matrix_id*/) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

std::shared_ptr<dumper::Field>
ContactMechanicsModel::createNodalFieldBool(const std::string & /*field_name*/,
                                            const std::string & /*group_name*/,
                                            bool /*padding_flag*/) {

  std::shared_ptr<dumper::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field>
ContactMechanicsModel::createNodalFieldReal(const std::string & field_name,
                                            const std::string & group_name,
                                            bool padding_flag) {

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["contact_force"] = this->internal_force;
  real_nodal_fields["blocked_dofs"] = this->blocked_dofs;
  real_nodal_fields["normals"] = this->normals;
  real_nodal_fields["tangents"] = this->tangents;
  real_nodal_fields["gaps"] = this->gaps;
  real_nodal_fields["previous_gaps"] = this->previous_gaps;
  real_nodal_fields["areas"] = this->nodal_area;

  if (padding_flag) {
    return this->mesh.createNodalField(real_nodal_fields[field_name],
                                       group_name, 3);
  } else {
    return this->mesh.createNodalField(real_nodal_fields[field_name],
                                       group_name);
  }

  std::shared_ptr<dumper::Field> field;
  return field;
}

#else
/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field>
ContactMechanicsModel::createNodalFieldReal(const std::string & /*field_name*/,
                                            const std::string & /*group_name*/,
                                            bool /*padding_flag*/) {
  return nullptr;
}

#endif

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::dump(const std::string & dumper_name) {
  mesh.dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::dump(const std::string & dumper_name, UInt step) {
  mesh.dump(dumper_name, step);
}

/* ------------------------------------------------------------------------- */
void ContactMechanicsModel::dump(const std::string & dumper_name, Real time,
                                 UInt step) {
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::dump() { mesh.dump(); }

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::dump(UInt step) { mesh.dump(step); }

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::dump(Real time, UInt step) {
  mesh.dump(time, step);
}

/* -------------------------------------------------------------------------- */
UInt ContactMechanicsModel::getNbData(
    const Array<Element> & elements, const SynchronizationTag & /*tag*/) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;

  for (const Element & el : elements) {
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::packData(CommunicationBuffer & /*buffer*/,
                                     const Array<Element> & /*elements*/,
                                     const SynchronizationTag & /*tag*/) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::unpackData(CommunicationBuffer & /*buffer*/,
                                       const Array<Element> & /*elements*/,
                                       const SynchronizationTag & /*tag*/) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt ContactMechanicsModel::getNbData(
    const Array<UInt> & dofs, const SynchronizationTag & /*tag*/) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;

  AKANTU_DEBUG_OUT();
  return size * dofs.size();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::packData(CommunicationBuffer & /*buffer*/,
                                     const Array<UInt> & /*dofs*/,
                                     const SynchronizationTag & /*tag*/) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::unpackData(CommunicationBuffer & /*buffer*/,
                                       const Array<UInt> & /*dofs*/,
                                       const SynchronizationTag & /*tag*/) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
