/**
 * @file   element.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  Element helper class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_HH_
#define AKANTU_ELEMENT_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Element                                                                    */
/* -------------------------------------------------------------------------- */
class Element {
public:
  ElementType type;
  Idx element;
  GhostType ghost_type;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline constexpr ElementKind kind() const;

  inline constexpr bool operator==(const Element & elem) const {
    return std::tie(type, element, ghost_type) ==
           std::tie(elem.type, elem.element, elem.ghost_type);
  }

  inline constexpr bool operator!=(const Element & elem) const {
    return std::tie(type, element, ghost_type) !=
           std::tie(elem.type, elem.element, elem.ghost_type);
  }

  inline constexpr bool operator<(const Element & rhs) const;
};

#if __cplusplus < 201703L
namespace {
  const Element ElementNull{_not_defined, Idx(-1), _casper};
} // namespace
#else
  inline constexpr Element ElementNull{_not_defined, Idx(-1), _casper};
#endif

/* -------------------------------------------------------------------------- */
inline constexpr bool Element::operator<(const Element & rhs) const {
  return ((rhs == ElementNull) ||
          std::tie(ghost_type, type, element) <
              std::tie(rhs.ghost_type, rhs.type, rhs.element));
}

} // namespace akantu

namespace std {
inline string to_string(const akantu::Element & _this) {
  if (_this == akantu::ElementNull) {
    return "ElementNull";
  }

  string str = "Element [" + to_string(_this.type) + ", " +
               to_string(_this.element) + ", " + to_string(_this.ghost_type) +
               "]";
  return str;
}
} // namespace std

namespace akantu {

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const Element & _this) {
  stream << std::to_string(_this);
  return stream;
}
} // namespace akantu

#endif /* AKANTU_ELEMENT_HH_ */
