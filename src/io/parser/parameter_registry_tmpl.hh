/**
 * @file   parameter_registry_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 09 2012
 * @date last modification: Fri Mar 27 2015
 *
 * @brief implementation of the templated part of ParameterRegistry class and
 * the derivated ones
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_error.hh"
#include "parameter_registry.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <string>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PARAMETER_REGISTRY_TMPL_HH__
#define __AKANTU_PARAMETER_REGISTRY_TMPL_HH__

namespace akantu {

namespace debug {
  class ParameterException : public Exception {
  public:
    ParameterException(const std::string & name, const std::string & message)
        : Exception(message), name(name) {}
    const std::string & name;
  };

  class ParameterUnexistingException : public ParameterException {
  public:
    ParameterUnexistingException(const std::string & name,
                                 const ParameterRegistry & registery)
        : ParameterException(name, "Parameter " + name +
                                       " does not exists in this scope") {
      auto && params = registery.listParameters();
      this->_info =
          std::accumulate(params.begin(), params.end(),
                          this->_info + "\n Possible parameters are: ",
                          [](auto && str, auto && param) {
                            static auto first = true;
                            auto && ret = str + (first ? " " : ", ") + param;
                            first = false;
                            return ret;
                          });
    }
  };

  class ParameterAccessRightException : public ParameterException {
  public:
    ParameterAccessRightException(const std::string & name,
                                  const std::string & perm)
        : ParameterException(name, "Parameter " + name + " is not " + perm) {}
  };

  class ParameterWrongTypeException : public ParameterException {
  public:
    ParameterWrongTypeException(const std::string name,
                                const std::type_info & wrong_type,
                                const std::type_info & type)
        : ParameterException(name, "Parameter " + name +
                                       " type error, cannot convert " +
                                       debug::demangle(type.name()) + " to " +
                                       debug::demangle(wrong_type.name())) {}
  };
} // namespace debug
/* -------------------------------------------------------------------------- */
template <typename T>
const ParameterTyped<T> & Parameter::getParameterTyped() const {
  try {
    const ParameterTyped<T> & tmp =
        dynamic_cast<const ParameterTyped<T> &>(*this);
    return tmp;
  } catch (std::bad_cast &) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterWrongTypeException(name, typeid(T), this->type()));
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> ParameterTyped<T> & Parameter::getParameterTyped() {
  try {
    ParameterTyped<T> & tmp = dynamic_cast<ParameterTyped<T> &>(*this);
    return tmp;
  } catch (std::bad_cast &) {
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterWrongTypeException(name, typeid(T), this->type()));
  }
}

/* ------------------------------------------------------------------------ */
template <typename T, typename V> void Parameter::set(const V & value) {
  if (!(isWritable()))
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "writable"));
  ParameterTyped<T> & typed_param = getParameterTyped<T>();
  typed_param.setTyped(value);
}

/* ------------------------------------------------------------------------ */
inline void Parameter::setAuto(__attribute__((unused))
                               const ParserParameter & value) {
  if (!(isParsable()))
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "parsable"));
}

/* -------------------------------------------------------------------------- */
template <typename T> const T & Parameter::get() const {
  if (!(isReadable()))
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "readable"));
  const ParameterTyped<T> & typed_param = getParameterTyped<T>();
  return typed_param.getTyped();
}

/* -------------------------------------------------------------------------- */
template <typename T> T & Parameter::get() {
  ParameterTyped<T> & typed_param = getParameterTyped<T>();
  if (!(isReadable()) || !(this->isWritable()))
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "accessible"));
  return typed_param.getTyped();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline Parameter::operator T() const {
  return this->get<T>();
}

/* -------------------------------------------------------------------------- */
template <typename T>
ParameterTyped<T>::ParameterTyped(std::string name, std::string description,
                                  ParameterAccessType param_type, T & param)
    : Parameter(name, description, param_type), param(param) {}

/* -------------------------------------------------------------------------- */
template <typename T>
template <typename V>
void ParameterTyped<T>::setTyped(const V & value) {
  param = value;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParameterTyped<T>::setAuto(const ParserParameter & value) {
  Parameter::setAuto(value);
  param = static_cast<T>(value);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<std::string>::setAuto(const ParserParameter & value) {
  Parameter::setAuto(value);
  param = value.getValue();
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<Vector<Real>>::setAuto(const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  Vector<Real> tmp = in_param;
  for (UInt i = 0; i < param.size(); ++i) {
    param(i) = tmp(i);
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<Matrix<Real>>::setAuto(const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  Matrix<Real> tmp = in_param;
  for (UInt i = 0; i < param.rows(); ++i) {
    for (UInt j = 0; j < param.cols(); ++j) {
      param(i, j) = tmp(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> const T & ParameterTyped<T>::getTyped() const {
  return param;
}
/* -------------------------------------------------------------------------- */
template <typename T> T & ParameterTyped<T>::getTyped() { return param; }

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParameterTyped<T>::printself(std::ostream & stream) const {
  Parameter::printself(stream);
  stream << param << std::endl;
}

/* -------------------------------------------------------------------------- */
template <>
inline void ParameterTyped<bool>::printself(std::ostream & stream) const {
  Parameter::printself(stream);
  stream << std::boolalpha << param << std::endl;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void ParameterRegistry::registerParam(std::string name, T & variable,
                                      ParameterAccessType type,
                                      const std::string & description) {
  std::map<std::string, Parameter *>::iterator it = params.find(name);
  if (it != params.end())
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterException(
        name, "Parameter named " + name + " already registered."));
  ParameterTyped<T> * param =
      new ParameterTyped<T>(name, description, type, variable);
  params[name] = param;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void ParameterRegistry::registerParam(std::string name, T & variable,
                                      const T & default_value,
                                      ParameterAccessType type,
                                      const std::string & description) {
  variable = default_value;
  registerParam(name, variable, type, description);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename V>
void ParameterRegistry::setMixed(const std::string & name, const V & value) {
  std::map<std::string, Parameter *>::iterator it = params.find(name);
  if (it == params.end()) {
    if (consisder_sub) {
      for (SubRegisteries::iterator it = sub_registries.begin();
           it != sub_registries.end(); ++it) {
        it->second->setMixed<T>(name, value);
      }
    } else {
      AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
    }
  } else {
    Parameter & param = *(it->second);
    param.set<T>(value);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void ParameterRegistry::set(const std::string & name, const T & value) {
  this->template setMixed<T>(name, value);
}

/* -------------------------------------------------------------------------- */
template <typename T> T & ParameterRegistry::get(const std::string & name) {
  Parameters::iterator it = params.find(name);
  if (it == params.end()) {
    if (consisder_sub) {
      for (SubRegisteries::iterator it = sub_registries.begin();
           it != sub_registries.end(); ++it) {
        try {
          return it->second->get<T>(name);
        } catch (...) {
        }
      }
    }

    // nothing was found not even in sub registries
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name));
  }

  Parameter & param = *(it->second);
  return param.get<T>();
}

/* -------------------------------------------------------------------------- */
const Parameter & ParameterRegistry::get(const std::string & name) const {
  Parameters::const_iterator it = params.find(name);
  if (it == params.end()) {
    if (consisder_sub) {
      for (SubRegisteries::const_iterator it = sub_registries.begin();
           it != sub_registries.end(); ++it) {
        try {
          return it->second->get(name);
        } catch (...) {
        }
      }
    }

    // nothing was found not even in sub registries
    AKANTU_CUSTOM_EXCEPTION(debug::ParameterUnexistingException(name, *this));
  }

  Parameter & param = *(it->second);
  return param;
}

/* -------------------------------------------------------------------------- */
namespace {
  namespace details {
    template <class T, class R, class Enable = void> struct CastHelper {
      static R convert(const T &) { throw std::bad_cast(); }
    };

    template <class T, class R>
    struct CastHelper<T, R,
                      std::enable_if_t<std::is_convertible<T, R>::value>> {
      static R convert(const T & val) { return val; }
    };
  } // namespace details
} // namespace

template <typename T> inline ParameterTyped<T>::operator Real() const {
  if (not isReadable())
    AKANTU_CUSTOM_EXCEPTION(
        debug::ParameterAccessRightException(name, "accessible"));
  return details::CastHelper<T, Real>::convert(param);
}

} // namespace akantu

#endif /* __AKANTU_PARAMETER_REGISTRY_TMPL_HH__ */
