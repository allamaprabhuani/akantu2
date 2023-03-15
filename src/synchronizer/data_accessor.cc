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

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void DataAccessor<Element>::packUnpackNodalDataHelper(
    Array<T> & data, CommunicationBuffer & buffer,
    const Array<Element> & elements, const Mesh & mesh) {
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

/* ------------------------------------------------------------------------ */
template <typename T, bool pack_helper>
void DataAccessor<Element>::packUnpackElementalDataHelper(
    ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & element, bool per_quadrature_point_data,
    const FEEngine & fem) {
  for (const auto & el : element) {
    auto nb_quad_per_elem =
        per_quadrature_point_data
            ? fem.getNbIntegrationPoints(el.type, el.ghost_type)
            : 1;
    auto nb_component = data_to_pack(el.type, el.ghost_type).getNbComponent();
    auto && data = data_to_pack.get(el, nb_component * nb_quad_per_elem);

    if (pack_helper) {
      buffer << data;
    } else {
      buffer >> data;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void DataAccessor<Idx>::packUnpackDOFDataHelper(Array<T> & data,
                                                CommunicationBuffer & buffer,
                                                const Array<Idx> & dofs) {
  T * data_ptr = data.data();
  for (const auto & dof : dofs) {
    if (pack_helper) {
      buffer << data_ptr[dof];
    } else {
      buffer >> data_ptr[dof];
    }
  }
}

/* -------------------------------------------------------------------------- */
#define DECLARE_HELPERS(T)                                                     \
  template void DataAccessor<Element>::packUnpackNodalDataHelper<T, false>(    \
      Array<T> & data, CommunicationBuffer & buffer,                           \
      const Array<Element> & elements, const Mesh & mesh);                     \
  template void DataAccessor<Element>::packUnpackNodalDataHelper<T, true>(     \
      Array<T> & data, CommunicationBuffer & buffer,                           \
      const Array<Element> & elements, const Mesh & mesh);                     \
  template void                                                                \
  DataAccessor<Element>::packUnpackElementalDataHelper<T, false>(              \
      ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,     \
      const Array<Element> & element, bool per_quadrature_point_data,          \
      const FEEngine & fem);                                                   \
  template void DataAccessor<Element>::packUnpackElementalDataHelper<T, true>( \
      ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,     \
      const Array<Element> & element, bool per_quadrature_point_data,          \
      const FEEngine & fem);                                                   \
  template void DataAccessor<Idx>::packUnpackDOFDataHelper<T, true>(           \
      Array<T> & data, CommunicationBuffer & buffer, const Array<Idx> & dofs); \
  template void DataAccessor<Idx>::packUnpackDOFDataHelper<T, false>(          \
      Array<T> & data, CommunicationBuffer & buffer, const Array<Idx> & dofs)

/* -------------------------------------------------------------------------- */
DECLARE_HELPERS(Real);
DECLARE_HELPERS(Idx);
DECLARE_HELPERS(bool);
/* -------------------------------------------------------------------------- */

} // namespace akantu
