/**
 * @file   coupler_solid_contact_explicit.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Thu May 22 2019
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model
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

  this->registerFEEngineObject<MyFEEngineType>("CouplerSolidContact", mesh,
					       Model::spatial_dimension);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("coupler_solid_contact", id, true);
  this->mesh.addDumpMeshToDumper("coupler_solid_contact", mesh, Model::spatial_dimension,
				 _not_ghost, _ek_regular);
#endif

  this->registerDataAccessor(*this);

  solid = new SolidMechanicsModel(mesh, Model::spatial_dimension,
				  "solid_mechanics_model", 0, this->dof_manager);
  contact = new ContactMechanicsModel(mesh, Model::spatial_dimension,
				  "contact_mechanics_model", 0, this->dof_manager);
  
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
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
CouplerSolidContact::getDefaultSolverID(const AnalysisMethod & method) {

  switch (method) {
  case _explicit_contact: {
    return std::make_tuple("explicit_contact", _tsst_static);
  }
  case _implicit_contact: {
    return std::make_tuple("implicit_contact", _tsst_static);
  }
  case _explicit_dynamic_contact: {
    return std::make_tuple("explicit_dynamic_contact", _tsst_dynamic_lumped);
    break;
  }  
  default:
    return std::make_tuple("unkown", _tsst_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType CouplerSolidContact::getDefaultSolverType() const {
  return _tsst_dynamic_lumped;
}

  
/* -------------------------------------------------------------------------- */
ModelSolverOptions CouplerSolidContact::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case _tsst_dynamic_lumped: {
    options.non_linear_solver_type = _nls_lumped;
    options.integration_scheme_type["displacement"] = _ist_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case _tsst_dynamic: {
    options.non_linear_solver_type = _nls_lumped;
    options.integration_scheme_type["displacement"] = _ist_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case _tsst_static: {
    options.non_linear_solver_type = _nls_newton_raphson_contact;
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

  // computes the internal forces
  this->assembleInternalForces();
  
  auto & internal_force = solid->getInternalForce();
  auto & external_force = solid->getExternalForce();
 
  auto & contact_force  = contact->getInternalForce();
  auto & contact_map    = contact->getContactMap();

  switch (method) {
  case _explicit_dynamic_contact:
  case _explicit_contact: {
    for (auto & pair: contact_map) {
      auto & connectivity = pair.second.connectivity;
      for (auto node : connectivity) {
	for (auto s : arange(spatial_dimension)) 
	  external_force(node, s) = contact_force(node, s);
      }
    }
    break;
  }
  default:
    break;
  }
    
  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement",
                                           external_force, 1);
  this->getDOFManager().assembleToResidual("displacement",
                                           internal_force, 1);
  switch (method) {
  case _implicit_contact: {
    this->getDOFManager().assembleToResidual("displacement",
					     contact_force, 1);
    break;
  }
  default:
    break;
  }
}
  
/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_IN();

  auto & internal_force = solid->getInternalForce();
  auto & external_force = solid->getExternalForce();
  
  auto & contact_force  = contact->getInternalForce();
  auto & contact_map    = contact->getContactMap();

  
  switch (method) {
  case _explicit_dynamic_contact:  
  case _explicit_contact: {
    for (auto & pair: contact_map) {
      auto & connectivity = pair.second.connectivity;
      for (auto node : connectivity) {
	for (auto s : arange(spatial_dimension)) 
	  external_force(node, s) = contact_force(node, s);
      }
    }
    break;
  }
  default:
    break;
  }
   
  if ("external" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement",
                                             external_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  if ("internal" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement",
                                             internal_force, 1);
    switch (method) {
    case _implicit_contact: {
      this->getDOFManager().assembleToResidual("displacement",
					       contact_force, 1);
      break;
    }
    default:
      break;
    }
    
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_CUSTOM_EXCEPTION(
      debug::SolverCallbackResidualPartUnknown(residual_part));

  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */
void CouplerSolidContact::beforeSolveStep() {}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::afterSolveStep() {}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::predictor() {

  switch (method) {
  case _explicit_dynamic_contact: {
    Array<Real> displacement(0, Model::spatial_dimension);

    Array<Real> current_positions(0, Model::spatial_dimension);
    auto positions = mesh.getNodes();
    current_positions.copy(positions);
  
    auto us = this->getDOFManager().getDOFs("displacement");
    //const auto deltas = this->getDOFManager().getSolution("displacement");
    const auto blocked_dofs = this->getDOFManager().getBlockedDOFs("displacement");

    for (auto && tuple : zip(make_view(us),
			     make_view(blocked_dofs),
			     make_view(current_positions))) {
      auto & u           = std::get<0>(tuple);
      const auto & bld   = std::get<1>(tuple);
      auto & cp          = std::get<2>(tuple);

      if (not bld)
	cp += u;   
    }

    contact->setPositions(current_positions); 
    contact->search();
    break;
  }
  default:
    break;
  }


}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::corrector() {

  switch (method) {
  case _implicit_contact:  
  case _explicit_contact: {
    Array<Real> displacement(0, Model::spatial_dimension);

    Array<Real> current_positions(0, Model::spatial_dimension);
    auto positions = mesh.getNodes();
    current_positions.copy(positions);
  
    auto us = this->getDOFManager().getDOFs("displacement");
    const auto deltas = this->getDOFManager().getSolution("displacement");
    const auto blocked_dofs = this->getDOFManager().getBlockedDOFs("displacement");

    for (auto && tuple : zip(make_view(us),
			     deltas,
			     make_view(blocked_dofs),
			     make_view(current_positions))) {
      auto & u           = std::get<0>(tuple);
      const auto & delta = std::get<1>(tuple);
      const auto & bld   = std::get<2>(tuple);
      auto & cp          = std::get<3>(tuple);

      if (not bld)
	cp += u + delta;   
    }

    contact->setPositions(current_positions); 
    contact->search();

    break;
  }
  default:
    break;
  }
    
  /*auto & internal_force = solid->getInternalForce();
  auto & external_force = solid->getExternalForce();

  
  
  std::stringstream filename;
  filename << "out" << "-00" << step << ".csv";
 
  std::ofstream outfile(filename.str());
 
  outfile << "x,gap,residual" << std::endl;
  
  auto & contact_map    = contact->getContactMap();
  for (auto & pair: contact_map) {
      auto & connectivity = pair.second.connectivity;
      auto node = connectivity(0);
      if (pair.second.gap > 0) {
	outfile << positions(node, 0) << "," << pair.second.gap << ","
		<< external_force(node, 1) + internal_force(node, 1)  << std::endl;
      }
  }

  outfile.close();
  step++;*/
}  

/* -------------------------------------------------------------------------- */
MatrixType CouplerSolidContact::getMatrixType(const ID & matrix_id) {

  if (matrix_id == "K")
    return _symmetric;
  if (matrix_id == "M") {
    return _symmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMatrix(const ID & matrix_id) {

  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    solid->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M") {
    solid->assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");
  
  solid->assembleInternalForces();  
  contact->assembleInternalForces();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  solid->assembleStiffnessMatrix();

  switch (method) {
  case _implicit_contact: {
    contact->assembleStiffnessMatrix();
    break;
  }
  default:
    break;
  }
     
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMassLumped() {
  solid->assembleMassLumped();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMass() {
  solid->assembleMass();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMassLumped(GhostType ghost_type)  {
  solid->assembleMassLumped(ghost_type);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMass(GhostType ghost_type) {
  solid->assembleMass(ghost_type);
}
  
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
  
  if (field_name == "contact_force" or field_name == "normals" or
      field_name == "gaps" or field_name == "previous_gaps" or
      field_name == "areas" or  field_name == "tangents")
    field = contact->createNodalFieldReal(field_name, group_name, padding_flag);
  else
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

/* --------------------------------------------------------------------------
 */
void CouplerSolidContact::dump(const std::string & dumper_name) {
  solid->onDump();
  mesh.dump(dumper_name);
}

/* --------------------------------------------------------------------------
 */
void CouplerSolidContact::dump(const std::string & dumper_name, UInt step) {
  solid->onDump();
  mesh.dump(dumper_name, step);
}

/* -------------------------------------------------------------------------
 */
void CouplerSolidContact::dump(const std::string & dumper_name, Real time,
                               UInt step) {
  solid->onDump();
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::dump() {
  solid->onDump();
  mesh.dump();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::dump(UInt step) {
  solid->onDump();
  mesh.dump(step);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::dump(Real time, UInt step) {
  solid->onDump();
  mesh.dump(time, step);
}
  
  
} // namespace akantu
