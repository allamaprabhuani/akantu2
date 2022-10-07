/**
 * @file   element_group.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Wed Dec 09 2020
 *
 * @brief  Stores information relevent to the notion of domain boundary and
 * surfaces.
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
#include "element_group.hh"
#include "aka_csr.hh"
#include "dumpable.hh"
#include "dumpable_inline_impl.hh"
#include "group_manager.hh"
#include "group_manager_inline_impl.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#if defined(AKANTU_COHESIVE_ELEMENT)
#include "cohesive_element_inserter.hh"
#endif
#include <algorithm>
#include <iterator>
#include <sstream>
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper_paraview.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
ElementGroup::ElementGroup(const std::string & group_name, const Mesh & mesh,
                           NodeGroup & node_group, Int dimension,
                           const std::string & id)
    : mesh(mesh), name(group_name), elements("elements", id),
      node_group(node_group), dimension(dimension) {
  AKANTU_DEBUG_IN();

  this->registerDumper<DumperParaview>("paraview_" + group_name, group_name,
                                       true);
  this->addDumpFilteredMesh(mesh, elements, node_group.getNodes(),
                            _all_dimensions);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementGroup::ElementGroup(const ElementGroup & /*other*/) = default;

/* -------------------------------------------------------------------------- */
void ElementGroup::clear() { elements.free(); }

/* -------------------------------------------------------------------------- */
void ElementGroup::clear(ElementType type, GhostType ghost_type) {
  if (elements.exists(type, ghost_type)) {
    elements(type, ghost_type).clear();
  }
}

/* -------------------------------------------------------------------------- */
bool ElementGroup::empty() const { return elements.empty(); }

/* -------------------------------------------------------------------------- */
void ElementGroup::append(const ElementGroup & other_group) {
  AKANTU_DEBUG_IN();

  node_group.append(other_group.node_group);

  /// loop on all element types in all dimensions
  for (auto ghost_type : ghost_types) {
    for (auto type : other_group.elementTypes(_ghost_type = ghost_type,
                     _element_kind = _ek_not_defined)) {
      const auto & other_elem_list = other_group.elements(type, ghost_type);
      auto nb_other_elem = other_elem_list.size();

      auto & elem_list = elements.alloc(0, 1, type, ghost_type);

      auto nb_elem = elem_list.size();
      /// append new elements to current list
      elem_list.resize(nb_elem + nb_other_elem);
      std::copy(other_elem_list.begin(), other_elem_list.end(),
                elem_list.begin() + nb_elem);
    }
  }

  this->optimize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT) {
    ;
  }

  stream << space << "ElementGroup [" << std::endl;
  stream << space << " + name: " << name << std::endl;
  stream << space << " + dimension: " << dimension << std::endl;
  elements.printself(stream, indent + 1);
  node_group.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ElementGroup::optimize() {
  // increasing the locality of data when iterating on the element of a group
  for (auto ghost_type : ghost_types) {
    for (auto type : elements.elementTypes(_ghost_type = ghost_type)) {
      auto & els = elements(type, ghost_type);
      std::sort(els.begin(), els.end());

      auto end = std::unique(els.begin(), els.end());
      els.resize(end - els.begin());
    }
  }

  node_group.optimize();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::fillFromNodeGroup() {
  CSR<Element> node_to_elem;
  MeshUtils::buildNode2Elements(this->mesh, node_to_elem, this->dimension);

  std::set<Element> seen;

  for (const auto & node : this->node_group) {
    for (const auto & element : node_to_elem.getRow(node)) {
      if (this->dimension != _all_dimensions &&
          this->dimension != Mesh::getSpatialDimension(element.type)) {
        continue;
      }

      if (seen.find(element) != seen.end()) {
        continue;
      }

      auto nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
      auto conn = mesh.getConnectivity(element);

      Int count = 0;
      for (auto n : conn) {
        count += (this->node_group.getNodes().find(n) != Idx(-1) ? 1 : 0);
      }

      if (count == nb_nodes_per_element) {
        this->add(element);
      }

      seen.insert(element);
    }
  }

  this->optimize();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::addDimension(Int dimension) {
  this->dimension = std::max(dimension, this->dimension);
}

/* -------------------------------------------------------------------------- */
void ElementGroup::onNodesAdded(const Array<Idx> & /*new_nodes*/,
                                const NewNodesEvent & event [[gnu::unused]]) {
#if defined(AKANTU_COHESIVE_ELEMENT)
  if (aka::is_of_type<CohesiveNewNodesEvent>(event)) {
    // nodes might have changed in the connectivity
    node_group.clear();
    const auto & mesh_to_mesh_facet =
        mesh.getData<Element>("mesh_to_mesh_facet");

    for (auto ghost_type : ghost_types) {
      for (auto type : elements.elementTypes(_ghost_type = ghost_type)) {
        auto & els = elements(type, ghost_type);

        if (not mesh_to_mesh_facet.exists(type, ghost_type)) {
          continue;
        }
        const auto & mesh_to_mesh_facet_type =
            mesh_to_mesh_facet(type, ghost_type);

        auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
        auto && conn_it = make_view(mesh.getConnectivity(type, ghost_type),
                                    nb_nodes_per_element)
                              .begin();

        auto && mesh_facet_conn_it =
            make_view(mesh.getMeshFacets().getConnectivity(type, ghost_type),
                      nb_nodes_per_element)
                .begin();

        for (auto element : els) {
          auto && mesh_facet_conn =
              mesh_facet_conn_it[mesh_to_mesh_facet_type(element).element];
          auto && conn = conn_it[element];
          conn = mesh_facet_conn;
          for (auto && n : conn) {
            node_group.add(n, false);
          }
        }
      }
    }
    node_group.optimize();
  }
#endif
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
