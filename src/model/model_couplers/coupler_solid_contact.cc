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
#include "dumpable_inline_impl.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

CouplerSolidContact::CouplerSolidContact(Mesh & mesh, UInt dim, const ID & id,
                                         std::shared_ptr<DOFManager> dof_manager,
                                         const ModelType model_type)
  : Model(mesh, model_type, dof_manager, dim, id) {

  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<MyFEEngineType>(
      "CouplerSolidContact", mesh, Model::spatial_dimension);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("coupler_solid_contact", id, true);
  this->mesh.addDumpMeshToDumper("coupler_solid_contact", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_regular);
#endif

  this->registerDataAccessor(*this);

  
  solid = new SolidMechanicsModel(mesh, Model::spatial_dimension,
				  "solid_mechanics_model",
				  0,
				  this->dof_manager);
  contact = new ContactMechanicsModel(mesh, Model::spatial_dimension,
				      "contact_mechanics_model",
				      0, 
				      this->dof_manager);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
CouplerSolidContact::~CouplerSolidContact() {}

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
void CouplerSolidContact::initSolver(TimeStepSolverType, NonLinearSolverType) {
  DOFManager & dof_manager = this->getDOFManager();

  /*this->allocNodalField(this->displacement, spatial_dimension, "displacement");
  this->allocNodalField(this->displacement_increment, spatial_dimension,
                       "displacement_increment");
  this->allocNodalField(this->external_force, spatial_dimension,
                        "external_force");
  if (not dof_manager.hasDOFs("displacement")) {
    dof_manager.registerDOFs("displacement", *this->displacement, _dst_nodal);
    }*/
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
    const TimeStepSolverType & type) const {
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

  solid->assembleInternalForces();
  contact->assembleInternalForces();
  
  auto & contact_force  = contact->getInternalForce();
  auto & external_force = solid->getExternalForce();
  auto & internal_force = solid->getInternalForce();
  
  /*auto & blocked_dofs   = solid->getBlockedDOFs();

  for (auto && values : zip(make_view(external_force),
			    make_view(contact_force),
                            make_view(blocked_dofs))) {
    auto & f_ext = std::get<0>(values);
    auto & f_con = std::get<1>(values);
    auto & blocked = std::get<2>(values);

    if (!blocked)
      f_ext = f_con;
  }*/

  
  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement",
                                           external_force, 1);
  this->getDOFManager().assembleToResidual("displacement",
                                           internal_force, 1);
  this->getDOFManager().assembleToResidual("displacement",
                                           contact_force, 1);   
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::beforeSolveStep() {
  contact->search();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::afterSolveStep() {}

/* -------------------------------------------------------------------------- */
MatrixType CouplerSolidContact::getMatrixType(const ID & matrix_id) {

  if (matrix_id == "K")
    return _symmetric;

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMatrix(const ID & matrix_id) {

  if (matrix_id == "K") {
    solid->assembleStiffnessMatrix();
    contact->assembleStiffnessMatrix();
  }

}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleLumpedMatrix(const ID & /*matrix_id*/) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::coupleExternalForces() {

  auto & contact_force = contact->getInternalForce();
  auto & external_force = solid->getExternalForce();
  auto & blocked_dofs = solid->getBlockedDOFs();

  for (auto && values : zip(make_view(external_force), make_view(contact_force),
                            make_view(blocked_dofs))) {
    auto & f_ext = std::get<0>(values);
    auto & f_con = std::get<1>(values);
    auto & blocked = std::get<2>(values);

    if (!blocked)
      f_ext = f_con;
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::coupleStiffnessMatrices() {
  auto & contact_stiffness =
      const_cast<SparseMatrix &>(contact->getDOFManager().getMatrix("K"));
  auto & solid_stiffness =
      const_cast<SparseMatrix &>(solid->getDOFManager().getMatrix("K"));

  solid_stiffness.add(contact_stiffness);
}

/* -------------------------------------------------------------------------- */
// void CouplerSolidContact::solve() {

// search for contacts
// contact.search();
// can be handled by beforeSolveStep()

/// assemble contact forces
/// contact.assembleInternalForces();
/// assemble contact stiffness matrix
/// contact.assembleStiffnessMatrix();
/// can be handled by solveStep but no need to solve it

/// couple the external forces
// this->coupleExternalForces();

// assemble the internal forces solid mechanics model
// assemble the stiffness forces
// couple the contact mechanics model stiffness to solid mechanics
// this->coupleStiffnessMatrices();

// and solve the solid mechanics model with new stifffness anf
// residual
/// all the above steps for solid mehcanics should be handled by the
// solveStep method with solvercall back to see that stifffness
// matrices are coupled

// ContactSolver callback(solid, contact);

// contact.solveStep();
// solid.solveStep(callback);
//}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
dumper::Field * CouplerSolidContact::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, const UInt & spatial_dimension,
    const ElementKind & kind) {

  dumper::Field * field = nullptr;

  field = solid->createElementalField(field_name, group_name, padding_flag,
				      spatial_dimension, kind);
  
  return field;
}


  
/* -------------------------------------------------------------------------- */
dumper::Field *
CouplerSolidContact::createNodalFieldReal(const std::string & field_name,
                                          const std::string & group_name,
                                          bool padding_flag) {

  dumper::Field * field = nullptr;
  
  field = solid->createNodalFieldReal(field_name, group_name, padding_flag);

  return field;
}

/* -------------------------------------------------------------------------- */  
dumper::Field * CouplerSolidContact::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  dumper::Field * field = nullptr;
  field = solid->createNodalFieldBool(field_name, group_name, padding_flag);
  return field;
}
  
#else

/* -------------------------------------------------------------------------- */
dumper::Field * CouplerSolidContact::createElementalField(const std::string &,
                                                          const std::string &,
                                                          bool, const UInt &,
                                                          const ElementKind &) {
  return nullptr;
}
/* ----------------------------------------------------------------------- */
dumper::Field * CouplerSolidContact::createNodalFieldReal(const std::string &,
                                                          const std::string &,
                                                          bool) {
  return nullptr;
}

/*-------------------------------------------------------------------*/
dumper::Field * CouplerSolidContact::createNodalFieldBool(const std::string &,
                                                          const std::string &,
                                                          bool) {
  return nullptr;
}
  
#endif

  
} // namespace akantu
