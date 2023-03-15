/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <memory>
#include <utility>

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Functions ArrayBase                                                       */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <> Idx Array<Real>::find(const Real & elem) const {
  AKANTU_DEBUG_IN();

  Real epsilon = std::numeric_limits<Real>::epsilon();
  auto it = std::find_if(begin(), end(), [&elem, &epsilon](auto && a) {
    return std::abs(a - elem) <= epsilon;
  });

  AKANTU_DEBUG_OUT();
  return (it != end()) ? end() - it : -1;
}

/* -------------------------------------------------------------------------- */
template <>
auto Array<ElementType>::operator*=(const ElementType & /*alpha*/) -> Array & {
  AKANTU_TO_IMPLEMENT();
}

template <>
auto Array<ElementType>::operator-=(const Array & /*vect*/) -> Array & {
  AKANTU_TO_IMPLEMENT();
}

template <>
auto Array<ElementType>::operator+=(const Array & /*vect*/) -> Array & {
  AKANTU_TO_IMPLEMENT();
}

template <> auto Array<char>::operator*=(const char & /*alpha*/) -> Array & {
  AKANTU_TO_IMPLEMENT();
}

template <> auto Array<char>::operator-=(const Array & /*vect*/) -> Array & {
  AKANTU_TO_IMPLEMENT();
}

template <> auto Array<char>::operator+=(const Array & /*vect*/) -> Array & {
  AKANTU_TO_IMPLEMENT();
}

} // namespace akantu
