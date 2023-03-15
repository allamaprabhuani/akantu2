/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_types.hh"
#include "element.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_QUADRATURE_POINT_H
#define AKANTU_QUADRATURE_POINT_H

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
class IntegrationPoint;
extern const IntegrationPoint IntegrationPointNull;
/* -------------------------------------------------------------------------- */

class IntegrationPoint : public Element {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IntegrationPoint(const Element & element = ElementNull, Int num_point = 0,
                   Int nb_quad_per_element = 0)
      : Element(element), num_point(num_point),
        global_num(element.element * nb_quad_per_element + num_point) {}

  constexpr IntegrationPoint(const IntegrationPoint & quad) = default;
  constexpr IntegrationPoint(IntegrationPoint && quad) = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline auto operator==(const IntegrationPoint & quad) const -> bool {
    return Element::operator==(quad) && this->num_point == quad.num_point;
  }

  inline auto operator!=(const IntegrationPoint & quad) const -> bool {
    return Element::operator!=(quad) || (num_point != quad.num_point) ||
           (global_num != quad.global_num);
  }

  auto operator<(const IntegrationPoint & rhs) const -> bool {
    bool res = Element::operator<(rhs) ||
               (Element::operator==(rhs) && this->num_point < rhs.num_point);
    return res;
  }

  inline auto operator=(const IntegrationPoint & q)
      -> IntegrationPoint & = default;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:
  /// number of quadrature point in the element
  Int num_point{0};
  /// global number of the quadrature point
  Int global_num{0};
};

} // namespace akantu

namespace std {
inline auto to_string(const akantu::IntegrationPoint & _this) -> string {
  if (static_cast<const akantu::Element &>(_this) == akantu::ElementNull) {
    return "IntegrationPoint[ElementNull]";
  }

  string str = "IntegrationPoint[ " +
               to_string(static_cast<const akantu::Element &>(_this)) + ", " +
               to_string(_this.num_point) + ": " + to_string(_this.global_num) +
               "]";
  return str;
}
} // namespace std

namespace akantu {

/// standard output stream operator
inline auto operator<<(std::ostream & stream, const IntegrationPoint & _this)
    -> std::ostream & {
  stream << std::to_string(_this);
  return stream;
}
} // namespace akantu

#endif /* AKANTU_QUADRATURE_POINT_H */
