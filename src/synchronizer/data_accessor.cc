/**
 * @file   data_accessor.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  data accessors constructor functions
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
  Int nb_nodes_per_element = 0;

  auto current_element_type = _not_defined;
  auto current_ghost_type = _casper;
  const_view_iterator<VectorProxy<Idx>> conn_it;

  auto data_it = make_view(data, nb_component).begin();

  for (const auto & el : elements) {
    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;
      nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
      conn_it = make_view(mesh.getConnectivity(el.type, el.ghost_type),
                          nb_nodes_per_element)
                    .begin();
    }

    auto && conn = conn_it[el.element];

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
  auto current_element_type = _not_defined;
  auto current_ghost_type = _casper;
  view_iterator<VectorProxy<T>> data_it;

  for (const auto & el : element) {
    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;
      auto nb_quad_per_elem =
          per_quadrature_point_data
              ? fem.getNbIntegrationPoints(el.type, el.ghost_type)
              : 1;
      auto && array = data_to_pack(el.type, el.ghost_type);
      auto nb_component = array.getNbComponent();
      data_it = make_view(array, nb_component * nb_quad_per_elem).begin();
    }

    auto && data = data_it[el.element];
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
