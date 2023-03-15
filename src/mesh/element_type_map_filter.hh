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

/* -------------------------------------------------------------------------- */
//#include "aka_array_filter.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_BY_ELEMENT_TYPE_FILTER_HH_
#define AKANTU_BY_ELEMENT_TYPE_FILTER_HH_
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
/* ElementTypeMapFilter */
/* -------------------------------------------------------------------------- */
template <class T, typename SupportType = ElementType>
class ElementTypeMapArrayFilter {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using array_type = ArrayFilter<T>;
  using value_type = typename array_type::value_type;

  using type_iterator =
      typename ElementTypeMapArray<Idx, SupportType>::type_iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ElementTypeMapArrayFilter(
      const ElementTypeMapArray<T, SupportType> & array,
      const ElementTypeMapArray<Idx, SupportType> & filter,
      const ElementTypeMap<Int, SupportType> & nb_data_per_elem)
      : array(array), filter(filter), nb_data_per_elem(nb_data_per_elem) {}

  ElementTypeMapArrayFilter(
      const ElementTypeMapArray<T, SupportType> & array,
      const ElementTypeMapArray<Idx, SupportType> & filter)
      : array(array), filter(filter) {}

  ~ElementTypeMapArrayFilter() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline const ArrayFilter<T>
  operator()(SupportType type, GhostType ghost_type = _not_ghost) const {
    if (filter.exists(type, ghost_type)) {
      Int nb_comp = 1;
      auto && array_v = array(type, ghost_type);
      if (nb_data_per_elem.exists(type, ghost_type)) {
        nb_comp = nb_data_per_elem(type, ghost_type) / array_v.getNbComponent();
      }
      return ArrayFilter<T>(array_v, filter(type, ghost_type), nb_comp);
    }
    return ArrayFilter<T>(empty_array, empty_filter, 1);
  };

  template <typename... Args>
  decltype(auto) elementTypes(Args &&... args) const {
    return filter.elementTypes(std::forward<decltype(args)>(args)...);
  }

  decltype(auto) getNbComponents(Int dim = _all_dimensions,
                                 GhostType ghost_type = _not_ghost,
                                 ElementKind kind = _ek_not_defined) const {
    return this->array.getNbComponents(dim, ghost_type, kind);
  };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  std::string getID() {
    return std::string("filtered:" + this->array().getID());
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:
  const ElementTypeMapArray<T, SupportType> & array;
  const ElementTypeMapArray<Idx, SupportType> & filter;
  ElementTypeMap<Int> nb_data_per_elem;

  /// Empty array to be able to return consistent filtered arrays
  Array<T> empty_array;
  Array<Idx> empty_filter;
};

} // namespace akantu

#endif /* AKANTU_BY_ELEMENT_TYPE_FILTER_HH_ */
