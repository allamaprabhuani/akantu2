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
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */


namespace akantu {

ContactMechanicsModel::ContactMechanicsModel( Mesh & mesh, UInt dim, const ID & id,
					      const MemoryID & memory_id,
					      const ModelType model_type)
  : Model(mesh, model_type, dim, id, memory_id) {

  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<MyFEEngineType>("ContactMechanicsModel", mesh,
					       Model::spatial_dimension);
  this->initDOFManager();

  this->detector = std::make_unique<ContactDetector>(this->mesh, id + ":contact_detector");

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
  
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::instantiateResolutions() {
  ParserSection model_section;
  bool is_empty;
  std::tie(model_section, is_empty) = this->getParserSection();

  if (not is_empty) {
    auto model_resolutions = model_section.getSubSections(ParserType::_contact_resolution);
    for (const auto & section : model_resolutions) {
      this->registerNewResolution(section);
    }
  }

  auto sub_sections = this->parser.getSubSections(ParserType::_contact_resolution);
  for (const auto & section : sub_sections) {
    this->registerNewResolution(section);
  }

  if (resolutions.empty())
    AKANTU_EXCEPTION("No contact resolutions where instantiated for the model"
                     << getID());
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
Resolution & ContactMechanicsModel::registerNewResolution(const ID & res_name,
							  const ID & res_type,
							  const ID & opt_param) {
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

  std::unique_ptr<Resolution> resolution = ResolutionFactory::getInstance().allocate(
      res_type, spatial_dimension, opt_param, *this, res_id);

  resolutions.push_back(std::move(resolution));

  return *(resolutions.back());
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::initResolutions() {
  AKANTU_DEBUG_ASSERT(resolutions.size() != 0, "No resolutions to initialize !");

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
std::tuple<ID, TimeStepSolverType>
ContactMechanicsModel::getDefaultSolverID(const AnalysisMethod & method) {

  /* @todo contact mechanics model doesnt really needs a solver have
     to find a way to fix this absurd part 
   */
  return std::make_tuple("contact", _tsst_static);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  /* ------------------------------------------------------------------------ */
  // computes the internal forces
  this->assembleInternalForces();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->contact_force, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the contact forces");

  // assemble the forces due to local stresses
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for (auto & resolution : resolutions) {
    resolution->assembleInternalForces(_not_ghost);
  }
  
  // assemble the stresses due to ghost elements
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for (auto & resolution : resolutions) {
    resolution->assembleInternalForces(_ghost);
  }

  AKANTU_DEBUG_OUT();

}

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
  // \TODO correct it for contact mechanics model, only one type of matrix
  if (matrix_id == "C")
    return _mt_not_defined;

  return _symmetric;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleMatrix(const ID & matrix_id) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::assembleLumpedMatrix(const ID & matrix_id) {
  AKANTU_TO_IMPLEMENT();
}
