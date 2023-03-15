/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_DUMPER_TYPE_TRAITS_HH_
#define AKANTU_DUMPER_TYPE_TRAITS_HH_
/* -------------------------------------------------------------------------- */
#include "element_type_map.hh"
#include "element_type_map_filter.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {
  /* ------------------------------------------------------------------------ */
  template <class data, class ret, class field> struct TypeTraits {
    //! the stored data (real, int, uint, ...)
    using data_type = data;
    //! the type returned by the operator *
    using return_type = ret;
    //! the field type (ElementTypeMap or ElementTypeMapFilter)
    using field_type = field;
    //! the type over which we iterate
    using it_type = typename field_type::value_type;
    //! the type of array (Array<T> or ArrayFilter<T>)
    using array_type = typename field_type::array_type;
    //! the iterator over the array
    using array_iterator = typename array_type::const_vector_iterator;
  };

  /* ------------------------------------------------------------------------ */

  // specialization for the case in which input and output types are the same
  template <class T, class ret, bool filtered>
  struct SingleType : public TypeTraits<T, ret, ElementTypeMapArray<T>> {};

  /* ------------------------------------------------------------------------ */

  // same as before but for filtered data
  template <class T, class ret>
  struct SingleType<T, ret, true>
      : public TypeTraits<T, ret, ElementTypeMapArrayFilter<T>> {};
  /* ------------------------------------------------------------------------ */

  // specialization for the case in which input and output types are different
  template <class it_type, class data_type, class ret, bool filtered>
  struct DualType
      : public TypeTraits<data_type, ret, ElementTypeMapArray<it_type>> {};

  /* ------------------------------------------------------------------------ */

  // same as before but for filtered data
  template <class it_type, class data_type, class ret>
  struct DualType<it_type, data_type, ret, true>
      : public TypeTraits<data_type, ret, ElementTypeMapArrayFilter<it_type>> {
  };
  /* ------------------------------------------------------------------------ */
} // namespace dumpers
} // namespace akantu

/* -------------------------------------------------------------------------- */

#endif /* AKANTU_DUMPER_TYPE_TRAITS_HH_ */
