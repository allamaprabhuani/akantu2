/**
 * @file   parser_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug  5 15:51:00 2013
 *
 * @brief  Implementation of the ParserParam conversions
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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <typename T>
inline ParserParameter::operator T() const {
  T t;
  std::stringstream sstr(value);
  sstr >> t;
  if(sstr.bad())
    AKANTU_EXCEPTION("No known conversion of a ParserParameter \""
                     << name << "\" to the type "
                     << typeid(T).name());
  return t;
}

/* -------------------------------------------------------------------------- */
template<>
inline ParserParameter::operator std::string() const {
  return value;
}

/* -------------------------------------------------------------------------- */
template<>
inline ParserParameter::operator Real() const {
  return Parser::parseReal(value, *parent_section);
}

/* --------------------------------------------------------- ----------------- */
template<>
inline ParserParameter::operator bool() const {
  bool b;
  std::stringstream sstr(value);
  sstr >> std::boolalpha >> b;
  if(sstr.fail()) {
    sstr.clear();
    sstr >> std::noboolalpha >> b;
  }
  return b;
}

/* --------------------------------------------------------- ----------------- */
template<>
inline ParserParameter::operator Vector<Real>() const {
  return Parser::parseVector(value, *parent_section);
}

/* --------------------------------------------------------- ----------------- */
template<>
inline ParserParameter::operator Matrix<Real>() const {
  return Parser::parseMatrix(value, *parent_section);
}

/* -------------------------------------------------------------------------- */
template<>
inline ParserParameter::operator RandomParameter<Real>() const {
  return Parser::parseRandomParameter(value, *parent_section);
}





__END_AKANTU__
