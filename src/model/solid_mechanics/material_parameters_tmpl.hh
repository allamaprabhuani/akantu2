/**
 * @file   material_parameters_tmpl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Aug  8 17:36:54 2012
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


__END_AKANTU__
#include <typeinfo>
__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
template<typename T>
const MaterialParamTyped<T> & MaterialParam::getMaterialParamTyped() const {
  try {
    const MaterialParamTyped<T> & tmp = dynamic_cast<const MaterialParamTyped<T> &>(*this);
    return tmp;
  } catch (...) {
    AKANTU_EXCEPTION("The parameter named " << name << " is of type "
		     << debug::demangle(typeid(T).name()) <<".");
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
MaterialParamTyped<T> & MaterialParam::getMaterialParamTyped() {
  try {
    MaterialParamTyped<T> & tmp = dynamic_cast<MaterialParamTyped<T> &>(*this);
    return tmp;
  } catch (...) {
    AKANTU_EXCEPTION("The parameter named " << name << " is of type "
		     << debug::demangle(typeid(T).name()) <<".");
  }
}

/* ------------------------------------------------------------------------ */
template<typename T>
void MaterialParam::set(T & value) {
  MaterialParamTyped<T> & typed_param = getMaterialParamTyped<T>();
  if(!(isWritable())) AKANTU_EXCEPTION("The parameter named " << name << " is not writable.");
  typed_param.setTyped(value);
}

/* -------------------------------------------------------------------------- */
template<typename T>
const T & MaterialParam::get() const {
  const MaterialParamTyped<T> & typed_param = getMaterialParamTyped<T>();
  if(!(isReadable())) AKANTU_EXCEPTION("The parameter named " << name << " is not readable.");
  return typed_param.getTyped();
}


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<typename T>
MaterialParamTyped<T>::MaterialParamTyped(std::string name, std::string description,
					  ParamAccessType param_type, T & param) :
  MaterialParam(name, description, param_type), param(param) {}

/* -------------------------------------------------------------------------- */
template<typename T>
void MaterialParamTyped<T>::setTyped(T & value) { param = value; }

/* -------------------------------------------------------------------------- */
template<typename T>
const T & MaterialParamTyped<T>::getTyped() const { return param;}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void MaterialParamTyped<T>::setParam(std::string value) {
  MaterialParam::setParam(value);
  std::stringstream sstr(value);
  sstr >> param;
}

/* -------------------------------------------------------------------------- */
template<>
inline void MaterialParamTyped<std::string>::setParam(std::string value) {
  MaterialParam::setParam(value);
  param = value;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void MaterialParamTyped<T>::printself(std::ostream & stream) const {
  MaterialParam::printself(stream);
  stream << param << std::endl;
}

template<>
inline void MaterialParamTyped<bool>::printself(std::ostream & stream) const {
  MaterialParam::printself(stream);
  stream << (param ? "true" : "false") << std::endl;
}


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<typename T>
void MaterialParameters::registerParam(std::string name, T & variable,
				       ParamAccessType type,
				       const std::string description) {
  std::map<std::string, MaterialParam *>::iterator it = params.find(name);
  if(it != params.end()) AKANTU_EXCEPTION("Parameter named " << name << " already registered.");
  MaterialParamTyped<T> * param = new MaterialParamTyped<T>(name, description, type,
							    variable);
  params[name] = param;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void MaterialParameters::registerParam(std::string name, T & variable,
				       T default_value,
				       ParamAccessType type,
				       const std::string description) {
  variable = default_value;
  registerParam(name, variable, type, description);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void MaterialParameters::set(std::string name, T value) {
  std::map<std::string, MaterialParam *>::iterator it = params.find(name);
  if(it == params.end()) AKANTU_EXCEPTION("No parameter named " << name << " in the material.");
  MaterialParam & param = *(it->second);
  param.set(value);
}

/* -------------------------------------------------------------------------- */
template<typename T>
const T & MaterialParameters::get(std::string name) const {
  std::map<std::string, MaterialParam *>::const_iterator it = params.find(name);
  if(it == params.end()) AKANTU_EXCEPTION("No parameter named " << name << " in the material.");
  const MaterialParam & param = *(it->second);
  return param.get<T>();
}
