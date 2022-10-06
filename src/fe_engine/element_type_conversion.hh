/**
 * @file   element_type_conversion.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  conversion between different types
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_TYPE_CONVERSION_HH_
#define AKANTU_ELEMENT_TYPE_CONVERSION_HH_

namespace akantu {

template <class InType, class OutType>
constexpr inline auto convertType(InType /*unused*/) {
  return OutType();
}

template <>
constexpr inline auto convertType<ElementType, ElementType>(ElementType type) {
  return type;
}

template <>
constexpr inline auto
convertType<ElementType, InterpolationType>(ElementType type) {
  return tuple_dispatch_with_default<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        return ElementClassProperty<type>::interpolation_type;
      },
      type, [&](auto && /*type*/) { return _itp_not_defined; });
}

} // namespace akantu

#endif /* AKANTU_ELEMENT_TYPE_CONVERSION_HH_ */
