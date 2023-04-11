/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include <regex>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T> inline ParserParameter::operator T() const {
  T t;
  std::stringstream sstr(value);
  sstr >> t;
  if (sstr.bad()) {
    AKANTU_EXCEPTION("No known conversion of a ParserParameter \""
                     << name << "\" to the type " << typeid(T).name());
  }
  return t;
}

#if !defined(DOXYGEN)
/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator const char *() const {
  return value.c_str();
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator Real() const {
  return Parser::parseReal(value, *parent_section);
}

/* --------------------------------------------------------- -----------------
 */
template <> inline ParserParameter::operator bool() const {
  bool b;
  std::stringstream sstr(value);
  sstr >> std::boolalpha >> b;
  if (sstr.fail()) {
    sstr.clear();
    sstr >> std::noboolalpha >> b;
  }
  return b;
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator std::vector<std::string>() const {
  std::vector<std::string> tmp;
  auto string =
      std::regex_replace(value, std::regex("[[:space:]]|\\[|\\]"), "");
  std::smatch sm;
  while (std::regex_search(string, sm, std::regex("[^,]+"))) {
    tmp.push_back(sm.str());
    string = sm.suffix();
  }
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator std::set<std::string>() const {
  std::set<std::string> tmp;
  auto string =
      std::regex_replace(value, std::regex("[[:space:]]|\\[|\\]"), "");
  std::smatch sm;
  while (std::regex_search(string, sm, std::regex("[^,]+"))) {
    tmp.emplace(sm.str());
    string = sm.suffix();
  }
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T, Int m, Int n, std::enable_if_t<n == 1> *>
inline ParserParameter::operator Matrix<T, m, n>() const {
  return Parser::parseVector(value, *parent_section);
}

// /* --------------------------------------------------------- ----------------
// */ template <Int m, Int n, std::enable_if_t<n == 1> * = nullptr> inline
// ParserParameter::operator Matrix<Int, m, n>() const {
//   Vector<Real> tmp = Parser::parseVector(value, *parent_section);
//   Vector<Int> tmp_int(tmp.size());
//   for (Int i = 0; i < tmp.size(); ++i) {
//     tmp_int(i) = Int(tmp(i));
//   }
//   return tmp_int;
// }

/* --------------------------------------------------------- ---------------- */
template <typename T, Int m, Int n, std::enable_if_t<n != 1> *>
inline ParserParameter::operator Matrix<T, m, n>() const {
  return Parser::parseMatrix(value, *parent_section);
}

/* -------------------------------------------------------------------------- */
template <> inline ParserParameter::operator RandomParameter<Real>() const {
  return Parser::parseRandomParameter(value, *parent_section);
}
#endif
} // namespace akantu
