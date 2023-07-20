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
#include "data_accessor.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DATA_ACCESSOR_TMPL_HH
#define AKANTU_DATA_ACCESSOR_TMPL_HH

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void DataAccessor<Element>::packUnpackNodalDataHelper(
    std::conditional_t<pack_helper, const Array<T>, Array<T>> & data,
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const Mesh & mesh) {
  Int nb_component = data.getNbComponent();
  auto data_it = make_view(data, nb_component).begin();

  for (const auto & el : elements) {
    auto && conn = mesh.getConnectivity(el);

    for (auto node : conn) {
      auto && data_vect = data_it[node];
      if (pack_helper) {
        buffer << data_vect;
      } else {
        buffer >> data_vect;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <
    typename T, bool pack_helper, class Func,
    std::enable_if_t<not std::is_base_of_v<FEEngine, std::decay_t<Func>>> *>
void DataAccessor<Element>::packUnpackElementalDataHelper(
    std::conditional_t<pack_helper, const ElementTypeMapArray<T>,
                       ElementTypeMapArray<T>> & data_to_pack,
    CommunicationBuffer & buffer, const Array<Element> & element,
    Func && data_per_element) {
  for (const auto & el : element) {
    auto nb_component = data_to_pack(el.type, el.ghost_type).getNbComponent();
    auto nb_data_per_elem = data_per_element(el) * nb_component;

    auto && data = data_to_pack.get(el, nb_data_per_elem);

    if (pack_helper) {
      buffer << data;
    } else {
      buffer >> data;
    }
  }
}
/* ------------------------------------------------------------------------ */
template <typename T, bool pack_helper>
void DataAccessor<Element>::packUnpackElementalDataHelper(
    std::conditional_t<pack_helper, const ElementTypeMapArray<T>,
                       ElementTypeMapArray<T>> & data_to_pack,
    CommunicationBuffer & buffer, const Array<Element> & element) {
  packUnpackElementalDataHelper<T, pack_helper>(
      data_to_pack, buffer, element, [](auto && /*el*/) { return 1; });
}

/* ------------------------------------------------------------------------ */
template <typename T, bool pack_helper>
void DataAccessor<Element>::packUnpackElementalDataHelper(
    std::conditional_t<pack_helper, const ElementTypeMapArray<T>,
                       ElementTypeMapArray<T>> & data_to_pack,
    CommunicationBuffer & buffer, const Array<Element> & element,
    const FEEngine & fem) {
  packUnpackElementalDataHelper<T, pack_helper>(
      data_to_pack, buffer, element, [&fem](auto && el) {
        return fem.getNbIntegrationPoints(el.type, el.ghost_type);
      });
}

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void DataAccessor<Idx>::packUnpackDOFDataHelper(
    std::conditional_t<pack_helper, const Array<T>, Array<T>> & data,
    CommunicationBuffer & buffer, const Array<Idx> & dofs) {
  auto data_ptr = make_view(data).begin();
  for (const auto & dof : dofs) {
    if constexpr (pack_helper) {
      buffer << data_ptr[dof];
    } else {
      buffer >> data_ptr[dof];
    }
  }
}

} // namespace akantu

#endif /* AKANTU_DATA_ACCESSOR_TMPL_HH */
