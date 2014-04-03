/**
 * @file   cppargparse_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Mar 31 21:15:31 2014
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
#include <stdexcept>

#ifndef __CPPARGPARSE_TMPL_HH__
#define __CPPARGPARSE_TMPL_HH__

namespace cppargparse {
/* -------------------------------------------------------------------------- */
/* Argument                                                                   */
/* -------------------------------------------------------------------------- */
struct ArgumentParser::_Argument : public Argument {
  _Argument() : Argument(),
		help(std::string()),
		nargs(1),
		type(_string),
		required(false),
		parsed(false),
		has_default(false),
		has_const(false),
		is_positional(false) {}
  virtual ~_Argument() {}

  void setValues(std::vector<std::string> & values) {
    for (std::vector<std::string>::iterator it = values.begin();
	 it != values.end(); ++it) {
      this->addValue(*it);
    }
  }

  virtual void addValue(std::string & value) = 0;
  virtual void setToDefault() = 0;
  virtual void setToConst() = 0;

  std::ostream & printDefault(std::ostream & stream) const {
    stream << std::boolalpha;
    if(has_default) {
      stream << " (default: ";
      this->_printDefault(stream);
      stream << ")";
    }
    if(has_const) {
      stream << " (const: ";
      this->_printConst(stream);
      stream << ")";
    }
    return stream;
  }

  virtual std::ostream & _printDefault(std::ostream & stream) const = 0;
  virtual std::ostream & _printConst(std::ostream & stream) const = 0;

  std::string help;
  int nargs;
  ArgumentType type;
  bool required;
  bool parsed;
  bool has_default;
  bool has_const;
  std::vector<std::string> keys;
  bool is_positional;
};

/* -------------------------------------------------------------------------- */
template<class T>
class ArgumentParser::ArgumentStorage : public ArgumentParser::_Argument {
public:
  ArgumentStorage() : _default(T()), _const(T()), values(std::vector<T>()) {}

  virtual void addValue(std::string & value) {
    std::stringstream sstr(value);
    T t;
    sstr >> t;
    values.push_back(t);
  }

  virtual void setToDefault() {
    values.clear();
    values.push_back(_default);
  }

  virtual void setToConst() {
    values.clear();
    values.push_back(_const);
  }

  void printself(std::ostream & stream) const {
    stream << this->name << " =";
    stream << std::boolalpha; // for boolean
    for(typename std::vector<T>::const_iterator vit = this->values.begin();
	vit != this->values.end(); ++vit) {
      stream << " " << *vit;
    }
  }

  virtual std::ostream & _printDefault(std::ostream & stream) const {
    stream << _default;
    return stream;
  }

  virtual std::ostream & _printConst(std::ostream & stream) const {
    stream << _const;
    return stream;
  }

  T _default;
  T _const;
  std::vector<T> values;
};

/* -------------------------------------------------------------------------- */
template<>
inline void ArgumentParser::ArgumentStorage<std::string>::addValue(std::string & value) {
  values.push_back(value);
}

template<class T>
struct is_vector {
  enum { value = false };
};


template<class T>
struct is_vector< std::vector<T> > {
  enum { value = true };
};


/* -------------------------------------------------------------------------- */
template<class T, bool is_vector = is_vector<T>::value>
struct cast_helper {
  static T cast(const ArgumentParser::Argument & arg) {
    const ArgumentParser::ArgumentStorage<T> & _arg =
      dynamic_cast<const ArgumentParser::ArgumentStorage<T> &>(arg);
    if(_arg.values.size() == 1) {
      return _arg.values[0];
    } else {
      throw std::length_error("Not enougth or too many argument where passed for the command line argument: " + arg.name);
    }
  }
};

template<class T>
struct cast_helper<T, true> {
  static T cast(const ArgumentParser::Argument & arg) {
    const ArgumentParser::ArgumentStorage<T> & _arg =
      dynamic_cast<const ArgumentParser::ArgumentStorage<T> &>(arg);
    return _arg.values;
  }
};

/* -------------------------------------------------------------------------- */
template<class T>
ArgumentParser::Argument::operator T() const{
  return cast_helper<T>::cast(*this);
}

template<class T>
void ArgumentParser::addArgument(const std::string & name_or_flag,
				 const std::string & help,
				 int nargs,
				 ArgumentType type,
				 T def) {
  _Argument & arg = _addArgument(name_or_flag, help, nargs, type);
  dynamic_cast<ArgumentStorage<T> &>(arg)._default = def;
  arg.has_default = true;
}

template<class T>
void ArgumentParser::addArgument(const std::string & name_or_flag,
				 const std::string & help,
				 int nargs,
				 ArgumentType type,
				 T def,
				 T cons) {
  _Argument & arg = _addArgument(name_or_flag, help, nargs, type);
  dynamic_cast<ArgumentStorage<T> &>(arg)._default = def;
  arg.has_default = true;
  dynamic_cast<ArgumentStorage<T> &>(arg)._const = cons;
  arg.has_const = true;
}


}

#endif /* __AKANTU_CPPARGPARSE_TMPL_HH__ */
