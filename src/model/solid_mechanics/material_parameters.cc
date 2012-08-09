/**
 * @file   material_parameters.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Aug  8 17:33:32 2012
 *
 * @brief  
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "material_parameters.hh"

__BEGIN_AKANTU__

MaterialParam::MaterialParam(): name(""), description(""), param_type(_pat_internal) {}

/* -------------------------------------------------------------------------- */
MaterialParam::MaterialParam(std::string name, std::string description,
			     ParamAccessType param_type) :
  name(name),
  description(description),
  param_type(param_type) {
  }

/* -------------------------------------------------------------------------- */
bool MaterialParam::isWritable() const { return param_type && _pat_writable; }

/* -------------------------------------------------------------------------- */
bool MaterialParam::isReadable() const { return param_type && _pat_readable; }

/* -------------------------------------------------------------------------- */
bool MaterialParam::isInternal() const { return param_type && _pat_internal; }

/* -------------------------------------------------------------------------- */
bool MaterialParam::isParsable() const { return param_type && _pat_parsable; }

/* -------------------------------------------------------------------------- */
void MaterialParam::printself(std::ostream & stream) const {
  std::stringstream sstr;
  sstr << name;
  UInt width = std::max(int(10 - sstr.str().length()), 0);
  sstr.width(width);

  if(description != "") {
    sstr << " [" << description << "]";
  }

  stream << sstr.str();
  width = std::max(int(50 - sstr.str().length()), 0);
  stream.width(width);

  stream << " : ";
}

/* -------------------------------------------------------------------------- */
void MaterialParam::setParam(std::string value) {
  if(!isParsable()) AKANTU_EXCEPTION("The parameter named " << name << " is not parsable.");
}


/* -------------------------------------------------------------------------- */
MaterialParameters::~MaterialParameters()  {
  std::map<std::string, MaterialParam *>::iterator it, end;
  for(it = params.begin(); it != params.end(); ++it)
    delete it->second;
}

/* -------------------------------------------------------------------------- */
void MaterialParameters::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  std::map<std::string, MaterialParam *>::const_iterator it, end;
  for(it = params.begin(); it != params.end(); ++it){
    stream << space << " + ";
    it->second->printself(stream);
  }
}

/* -------------------------------------------------------------------------- */
void MaterialParameters::setParam(std::string name, std::string value) {
  std::map<std::string, MaterialParam *>::iterator it = params.find(name);
  if(it == params.end()) AKANTU_EXCEPTION("No parameter named " << name << " registered.");

  MaterialParam & param = *(it->second);
  param.setParam(value);
}


__END_AKANTU__
