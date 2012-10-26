/**
 * @file   material_parameters.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Aug  8 11:30:15 2012
 *
 * @brief  class to handle the material parameters
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
#include "aka_common.hh"
#include <map>

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_PARAMETERS_HH__
#define __AKANTU_MATERIAL_PARAMETERS_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
enum ParamAccessType {
  _pat_internal   = 0x0001,
  _pat_writable   = 0x0010,
  _pat_readable   = 0x0100,
  _pat_modifiable = 0x0110, //_pat_readable | _pat_writable,
  _pat_parsable   = 0x1000,
  _pat_parsmod    = 0x1110
};

template<typename T> class MaterialParamTyped;
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
class MaterialParam {
public:
  MaterialParam();
  MaterialParam(std::string name, std::string description, ParamAccessType param_type);

  virtual ~MaterialParam() {};
  /* ------------------------------------------------------------------------ */
  bool isInternal() const;
  bool isWritable() const;
  bool isReadable() const;
  bool isParsable() const;

  /* ------------------------------------------------------------------------ */
  template<typename T> void set(T & value);
  template<typename T> const T & get() const;

  virtual void parseParam(std::string value);

  /* ------------------------------------------------------------------------ */
  virtual void printself(std::ostream & stream) const;

protected:
  template<typename T>
  const MaterialParamTyped<T> & getMaterialParamTyped() const;

  template<typename T>
  MaterialParamTyped<T> & getMaterialParamTyped();

private:
  std::string name;
  std::string description;
  ParamAccessType param_type;
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<typename T>
class MaterialParamTyped : public MaterialParam {
public:
  MaterialParamTyped(std::string name, std::string description,
		     ParamAccessType param_type, T & param);

  /* ------------------------------------------------------------------------ */
  void setTyped(T & value);
  const T & getTyped() const;

  void parseParam(std::string value);

  virtual void printself(std::ostream & stream) const;
private:
  T & param;
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
class MaterialParameters {
public:
  ~MaterialParameters();

  template<typename T>
  void registerParam(std::string name, T & variable,
		     ParamAccessType type,
		     const std::string description = "");

  template<typename T>
  void registerParam(std::string name, T & variable, T default_value,
		     ParamAccessType type,
		     const std::string description = "");

  /* ------------------------------------------------------------------------ */
  template<typename T> void set(std::string name, T value);
  template<typename T> const T & get(std::string name) const;

  void parseParam(std::string name, std::string value);

  /* ------------------------------------------------------------------------ */
  void printself(std::ostream & stream, int indent) const;

private:
  std::map<std::string, MaterialParam *> params;
};

#include "material_parameters_tmpl.hh"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_PARAMETERS_HH__ */
