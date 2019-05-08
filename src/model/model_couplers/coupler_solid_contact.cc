/**
 * @file   coupler_solid_contact_explicit.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * 
 * @date creation: Thu Jan 17 2019
 * @date last modification: Thu Jan 17 2019
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model in explicit 
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
#include "coupler_solid_contact.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "dumpable_inline_impl.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

CouplerSolidContact::CouplerSolidContact(SolidMechanicsModel & solid,
					 ContactMechanicsModel & contact, UInt dim,
					 const ID & id,
					 const ModelType model_type)
  : Model(solid.getMesh(),model_type, dim, id), 
    solid(solid), contact(contact), displacement(nullptr),
    displacement_increment(nullptr), external_force(nullptr) {

  AKANTU_DEBUG_IN();

 this->registerFEEngineObject<MyFEEngineType>("CouplerSolidContact", solid.getMesh(),
					      Model::spatial_dimension);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("coupler_solid_contact", id, true);
  this->mesh.addDumpMeshToDumper("coupler_solid_contact", solid.getMesh(),
			       Model::spatial_dimension, _not_ghost,
			       _ek_regular);
#endif

  this->initDOFManager();
  this->registerDataAccessor(*this);

  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
CouplerSolidContact::~CouplerSolidContact() {
 
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::initFullImpl(const ModelOptions & options) {

  Model::initFullImpl(options);

  this->initBC(*this, *displacement, *displacement_increment, *external_force);
}
  
/* -------------------------------------------------------------------------- */
void CouplerSolidContact::initModel() {
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
FEEngine & CouplerSolidContact::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(
      getFEEngineClassBoundary<MyFEEngineType>(name));
}

  
/* -------------------------------------------------------------------------- */
void CouplerSolidContact::initSolver(TimeStepSolverType time_step_solver_type,
							       NonLinearSolverType) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->displacement, spatial_dimension,
			"displacement");
  this->allocNodalField(this->displacement_increment, spatial_dimension,
                        "displacement_increment");
  this->allocNodalField(this->external_force, spatial_dimension,
                        "external_force");
  if (!dof_manager.hasDOFs("displacement")) {
    dof_manager.registerDOFs("displacement", *this->displacement, _dst_nodal);
  }
}
  
  
/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
CouplerSolidContact::getDefaultSolverID(const AnalysisMethod & method) {
  
  switch (method) {
  case _explicit_contact: {
    return std::make_tuple("explicit_contact", _tsst_dynamic);
  }
  case _implicit_contact: {
    return std::make_tuple("implicit_contact", _tsst_static);
  }  
  default:
    return std::make_tuple("unkown", _tsst_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions CouplerSolidContact::getDefaultSolverOptions(
   const TimeStepSolverType & type ) const {
  ModelSolverOptions options;
  
  switch (type) {
  case _tsst_dynamic: {
    options.non_linear_solver_type = _nls_newton_raphson;
    options.integration_scheme_type["displacement"] = _ist_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  case _tsst_static: {
    options.non_linear_solver_type = _nls_newton_raphson;
    options.integration_scheme_type["displacement"] = _ist_pseudo_time;
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
void CouplerSolidContact::assembleResidual() {
  solid.assembleInternalForces();
  contact.assembleInternalForces();

  this->coupleExternalForces();  
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::beforeSolveStep() {
  contact.search();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::afterSolveStep() {
  
}

  
/* -------------------------------------------------------------------------- */
MatrixType CouplerSolidContact::getMatrixType(const ID & matrix_id) {
 
  if (matrix_id == "K")
    return _symmetric;

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMatrix(const ID & matrix_id) {

  contact.assembleStiffnessMatrix();

  auto & solid_tss = solid.getTimeStepSolver();
  //solid.assembleMatrix(matrix_id);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleLumpedMatrix(const ID & matrix_id) {
  AKANTU_TO_IMPLEMENT();
}
								     
 
/* -------------------------------------------------------------------------- */
void CouplerSolidContact::coupleExternalForces() {

  auto & contact_force  = contact.getInternalForce();
  auto & external_force = solid.getExternalForce();
  auto & blocked_dofs   = solid.getBlockedDOFs();
  
  for (auto && values: zip(make_view(external_force),
			   make_view(contact_force),
			   make_view(blocked_dofs)) ) {
    auto & f_ext   = std::get<0>(values);
    auto & f_con   = std::get<1>(values);
    auto & blocked = std::get<2>(values);

    if (!blocked) 
      f_ext = f_con;
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::coupleStiffnessMatrices() {
  auto & contact_stiffness =
    const_cast<SparseMatrix &>(contact.getDOFManager().getMatrix("K"));
  auto & solid_stiffness =
    const_cast<SparseMatrix &>(solid.getDOFManager().getMatrix("K"));
  
  solid_stiffness.add(contact_stiffness);
}

/* -------------------------------------------------------------------------- */
//void CouplerSolidContact::solve() {

  // search for contacts
  // contact.search();
  // can be handled by beforeSolveStep()

  /// assemble contact forces
  /// contact.assembleInternalForces();
  /// assemble contact stiffness matrix
  /// contact.assembleStiffnessMatrix();
  /// can be handled by solveStep but no need to solve it
  
  /// couple the external forces
  //this->coupleExternalForces();

  // assemble the internal forces solid mechanics model
  // assemble the stiffness forces
  // couple the contact mechanics model stiffness to solid mechanics
  //this->coupleStiffnessMatrices();

  // and solve the solid mechanics model with new stifffness anf
  // residual
  /// all the above steps for solid mehcanics should be handled by the
  // solveStep method with solvercall back to see that stifffness
  // matrices are coupled
    

  //ContactSolver callback(solid, contact);
  
  //contact.solveStep();
  //solid.solveStep(callback);
  //}

/* -------------------------------------------------------------------------- */
UInt CouplerSolidContact::getNbData(const Array<Element> & elements,
                                    const SynchronizationTag & tag) const {
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
void CouplerSolidContact::packData(CommunicationBuffer & buffer,
                                   const Array<Element> & elements,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::unpackData(CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt CouplerSolidContact::getNbData(const Array<UInt> & dofs,
                                    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  AKANTU_DEBUG_OUT();
  return size * dofs.size();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::packData(CommunicationBuffer & buffer,
                                   const Array<UInt> & dofs,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::unpackData(CommunicationBuffer & buffer,
                                     const Array<UInt> & dofs,
                                     const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
 
dumper::Field *
CouplerSolidContact::createNodalFieldBool(const std::string & field_name,
					    const std::string & group_name,
					    __attribute__((unused)) bool padding_flag) {

  dumper::Field * field = nullptr;
  return field;
}

/* -------------------------------------------------------------------------- */  
dumper::Field *
CouplerSolidContact::createNodalFieldReal(const std::string & field_name,
					  const std::string & group_name,
					  bool padding_flag) {

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["displacement"]   = this->displacement;
  real_nodal_fields["external_force"] = this->external_force;
    
  dumper::Field * field = nullptr;
  if (padding_flag)
    field = this->mesh.createNodalField(real_nodal_fields[field_name],
                                        group_name, 3);
  else
    field =
        this->mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}
  
#else
/* -------------------------------------------------------------------------- */
dumper::Field * CouplerSolidContact::createNodalFieldReal(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}

#endif  
  
  
} // akantu
  
