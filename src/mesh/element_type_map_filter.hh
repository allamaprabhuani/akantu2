/**
 * @file   element_type_map_filter.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Filtered version based on a an akantu::ElementGroup of a
 * akantu::ElementTypeMap
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
#include "aka_array_filter.hh"
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
      typename ElementTypeMapArray<UInt, SupportType>::type_iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ElementTypeMapArrayFilter(
      const ElementTypeMapArray<T, SupportType> & array,
      const ElementTypeMapArray<UInt, SupportType> & filter,
      const ElementTypeMap<UInt, SupportType> & nb_data_per_elem)
      : array(array), filter(filter), nb_data_per_elem(nb_data_per_elem) {}

  ElementTypeMapArrayFilter(
      const ElementTypeMapArray<T, SupportType> & array,
      const ElementTypeMapArray<UInt, SupportType> & filter)
      : array(array), filter(filter) {}

  ~ElementTypeMapArrayFilter() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  inline ArrayFilter<T> operator()(const SupportType & type,
                                   GhostType ghost_type = _not_ghost) const {
    if (filter.exists(type, ghost_type)) {
      if (nb_data_per_elem.exists(type, ghost_type)) {
        return ArrayFilter<T>(array(type, ghost_type), filter(type, ghost_type),
                              nb_data_per_elem(type, ghost_type) /
                                  array(type, ghost_type).getNbComponent());
      }
      return ArrayFilter<T>(array(type, ghost_type), filter(type, ghost_type),
                            1);
    }
    return ArrayFilter<T>(empty_array, empty_filter, 1);
  };

  template <typename... Args>
  decltype(auto) elementTypes(Args &&... args) const {
    return filter.elementTypes(std::forward<decltype(args)>(args)...);
  }

  decltype(auto) getNbComponents(UInt dim = _all_dimensions,
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
  const ElementTypeMapArray<UInt, SupportType> & filter;
  ElementTypeMap<UInt> nb_data_per_elem;

  /// Empty array to be able to return consistent filtered arrays
  Array<T> empty_array;
  Array<UInt> empty_filter;
};

} // namespace akantu

#endif /* AKANTU_BY_ELEMENT_TYPE_FILTER_HH_ */
