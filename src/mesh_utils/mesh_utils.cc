/**
 * @file   mesh_utils.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Aug 20 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  All mesh utils necessary for various tasks
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "mesh_utils.hh"
#include "element_synchronizer.hh"
#include "fe_engine.hh"
#include "mesh_accessor.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <limits>
#include <numeric>
#include <queue>
#include <set>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNode2Elements(const Mesh & mesh,
                                   CSR<Element> & node_to_elem,
                                   UInt spatial_dimension,
                                   ElementKind el_kind) {
  AKANTU_DEBUG_IN();
  if (spatial_dimension == _all_dimensions) {
    spatial_dimension = mesh.getSpatialDimension();
  }

  /// count number of occurrence of each node
  UInt nb_nodes = mesh.getNbNodes();

  /// array for the node-element list
  node_to_elem.resizeRows(nb_nodes);
  node_to_elem.clearRows();

  for_each_element(
      mesh,
      [&](auto && element) {
        Vector<UInt> conn = mesh.getConnectivity(element);
        std::set<UInt> unique_nodes;
        for (auto && node : conn) {
          auto ret = unique_nodes.emplace(node);
          if (ret.second) {
            ++node_to_elem.rowOffset(node);
          }
        }
      },
      _spatial_dimension = spatial_dimension, _element_kind = el_kind);

  node_to_elem.countToCSR();
  node_to_elem.resizeCols();

  /// rearrange element to get the node-element list
  // Element e;
  node_to_elem.beginInsertions();

  for_each_element(
      mesh,
      [&](auto && element) {
        Vector<UInt> conn = mesh.getConnectivity(element);
        std::set<UInt> unique_nodes;
        for (auto && node : conn) {
          auto ret = unique_nodes.emplace(node);
          if (ret.second) {
            node_to_elem.insertInRow(node, element);
          }
        }
      },
      _spatial_dimension = spatial_dimension, _element_kind = el_kind);

  node_to_elem.endInsertions();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNode2ElementsElementTypeMap(const Mesh & mesh,
                                                 CSR<UInt> & node_to_elem,
                                                 ElementType type,
                                                 GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_elements = mesh.getConnectivity(type, ghost_type).size();

  UInt * conn_val = mesh.getConnectivity(type, ghost_type).storage();

  /// array for the node-element list
  node_to_elem.resizeRows(nb_nodes);
  node_to_elem.clearRows();

  /// count number of occurrence of each node
  for (UInt el = 0; el < nb_elements; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      ++node_to_elem.rowOffset(conn_val[el_offset + n]);
    }
  }

  /// convert the occurrence array in a csr one
  node_to_elem.countToCSR();

  node_to_elem.resizeCols();
  node_to_elem.beginInsertions();

  /// save the element index in the node-element list
  for (UInt el = 0; el < nb_elements; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      node_to_elem.insertInRow(conn_val[el_offset + n], el);
    }
  }

  node_to_elem.endInsertions();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildFacets(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (auto ghost_type : ghost_types) {
    for (const auto & type :
         mesh.elementTypes(spatial_dimension - 1, ghost_type)) {
      mesh.getConnectivity(type, ghost_type).resize(0);
      // \todo inform the mesh event handler
    }
  }

  buildFacetsDimension(mesh, mesh, true, spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildAllFacets(const Mesh & mesh, Mesh & mesh_facets,
                               UInt to_dimension) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  buildAllFacets(mesh, mesh_facets, spatial_dimension, to_dimension);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void MeshUtils::buildAllFacets(const Mesh & mesh, Mesh & mesh_facets,
                               UInt from_dimension, UInt to_dimension) {
  AKANTU_DEBUG_IN();

  to_dimension = std::max(to_dimension, UInt(0));

  AKANTU_DEBUG_ASSERT(
      mesh_facets.isMeshFacets(),
      "The mesh_facets should be initialized with initMeshFacets");

  /// generate facets
  buildFacetsDimension(mesh, mesh_facets, false, from_dimension);

  /// sort facets and generate sub-facets
  for (UInt i = from_dimension - 1; i > to_dimension; --i) {
    buildFacetsDimension(mesh_facets, mesh_facets, false, i);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildFacetsDimension(const Mesh & mesh, Mesh & mesh_facets,
                                     bool boundary_only, UInt dimension) {
  AKANTU_DEBUG_IN();

  // save the current parent of mesh_facets and set it temporarly to mesh since
  // mesh is the one containing the elements for which mesh_facets has the
  // sub-elements
  // example: if the function is called with mesh = mesh_facets
  const Mesh * mesh_facets_parent = nullptr;
  try {
    mesh_facets_parent = &mesh_facets.getMeshParent();
  } catch (...) {
  }

  mesh_facets.defineMeshParent(mesh);
  MeshAccessor mesh_accessor(mesh_facets);

  UInt spatial_dimension = mesh.getSpatialDimension();

  const Array<Real> & mesh_facets_nodes = mesh_facets.getNodes();
  const auto mesh_facets_nodes_it = mesh_facets_nodes.begin(spatial_dimension);

  CSR<Element> node_to_elem;
  buildNode2Elements(mesh, node_to_elem, dimension);

  Array<UInt> counter;
  std::vector<Element> connected_elements;

  NewElementsEvent event(AKANTU_CURRENT_FUNCTION);

  // init the SubelementToElement data to improve performance
  for (auto && ghost_type : ghost_types) {
    for (auto && type : mesh.elementTypes(dimension, ghost_type)) {
      auto & subelement_to_element =
          mesh_accessor.getSubelementToElement(type, ghost_type);
      subelement_to_element.resize(mesh.getNbElement(type, ghost_type),
                                   ElementNull);

      auto facet_types = mesh.getAllFacetTypes(type);

      for (auto && ft : arange(facet_types.size())) {
        auto facet_type = facet_types(ft);
        mesh_accessor.getElementToSubelement(facet_type, ghost_type);
        mesh_accessor.getConnectivity(facet_type, ghost_type);
      }
    }
  }

  const ElementSynchronizer * synchronizer = nullptr;
  if (mesh.isDistributed()) {
    synchronizer = &(mesh.getElementSynchronizer());
  }

  Element current_element;
  for (auto && ghost_type : ghost_types) {
    GhostType facet_ghost_type = ghost_type;
    current_element.ghost_type = ghost_type;

    for (auto && type : mesh.elementTypes(dimension, ghost_type)) {
      auto facet_types = mesh.getAllFacetTypes(type);
      current_element.type = type;

      for (auto && ft : arange(facet_types.size())) {
        auto facet_type = facet_types(ft);
        auto nb_element = mesh.getNbElement(type, ghost_type);

        auto && element_to_subelement =
            &mesh_accessor.getElementToSubelementNC(facet_type, ghost_type);
        auto && connectivity_facets =
            &mesh_accessor.getConnectivity(facet_type, ghost_type);

        auto nb_nodes_per_facet = connectivity_facets->getNbComponent();
        // Vector<UInt> facet(nb_nodes_per_facet);

        for (UInt el = 0; el < nb_element; ++el) {
          current_element.element = el;

          auto && facets =
              mesh.getFacetConnectivity(current_element, ft).transpose();

          for (auto facet : facets) {
            // facet = facets(f);

            UInt first_node_nb_elements = node_to_elem.getNbCols(facet(0));
            counter.resize(first_node_nb_elements);
            counter.zero();

            // loop over the other nodes to search intersecting elements,
            // which are the elements that share another node with the
            // starting element after first_node
            for (auto && data : enumerate(node_to_elem.getRow(facet(0)))) {
              auto && local_el = std::get<0>(data);
              auto && first_node = std::get<1>(data);
              for (auto n : arange(1, nb_nodes_per_facet)) {
                auto && node_elements = node_to_elem.getRow(facet(n));
                counter(local_el) += std::count(
                    node_elements.begin(), node_elements.end(), first_node);
              }
            }

            // counting the number of elements connected to the facets and
            // taking the minimum element number, because the facet should
            // be inserted just once
            UInt nb_element_connected_to_facet = 0;
            Element minimum_el = ElementNull;
            connected_elements.clear();
            for (auto && data : enumerate(node_to_elem.getRow(facet(0)))) {

              if (not(counter(std::get<0>(data)) == nb_nodes_per_facet - 1)) {
                continue;
              }

              auto && real_el = std::get<1>(data);

              ++nb_element_connected_to_facet;
              minimum_el = std::min(minimum_el, real_el);
              connected_elements.push_back(real_el);
            }

            if (minimum_el != current_element) {
              continue;
            }

            bool full_ghost_facet = false;

            UInt n = 0;
            while (n < nb_nodes_per_facet and mesh.isPureGhostNode(facet(n))) {
              ++n;
            }
            if (n == nb_nodes_per_facet) {
              full_ghost_facet = true;
            }

            if (full_ghost_facet) {
              continue;
            }

            if (boundary_only and nb_element_connected_to_facet != 1) {
              continue;
            }

            std::vector<Element> elements;

            // build elements_on_facets: linearized_el must come first
            // in order to store the facet in the correct direction
            // and avoid to invert the sign in the normal computation
            elements.reserve(elements.size() + connected_elements.size());
            for (auto && connected_element : connected_elements) {
              elements.push_back(connected_element);
            }

            if (nb_element_connected_to_facet == 1) { /// boundary facet
              elements.push_back(ElementNull);
            } else if (nb_element_connected_to_facet == 2) { /// internal facet
              /// check if facet is in between ghost and normal
              /// elements: if it's the case, the facet is either
              /// ghost or not ghost. The criterion to decide this
              /// is arbitrary. It was chosen to check the processor
              /// id (prank) of the two neighboring elements. If
              /// prank of the ghost element is lower than prank of
              /// the normal one, the facet is not ghost, otherwise
              /// it's ghost
              GhostType gt[2] = {_not_ghost, _not_ghost};

              for (UInt el = 0; el < connected_elements.size(); ++el) {
                gt[el] = connected_elements[el].ghost_type;
              }

              if ((gt[0] == _not_ghost) xor (gt[1] == _not_ghost)) {
                UInt prank[2];
                for (UInt el = 0; el < 2; ++el) {
                  prank[el] = synchronizer->getRank(connected_elements[el]);
                }

                // ugly trick from Marco detected :P
                bool ghost_one = (gt[0] != _ghost);
                if (prank[ghost_one] > prank[!ghost_one]) {
                  facet_ghost_type = _not_ghost;
                } else {
                  facet_ghost_type = _ghost;
                }

                connectivity_facets = &mesh_accessor.getConnectivity(
                    facet_type, facet_ghost_type);
                element_to_subelement = &mesh_accessor.getElementToSubelementNC(
                    facet_type, facet_ghost_type);
              }
            }

            element_to_subelement->push_back(elements);
            connectivity_facets->push_back(facet);

            /// current facet index
            UInt current_facet = connectivity_facets->size() - 1;
            Element facet_element{facet_type, current_facet, facet_ghost_type};
            event.getList().push_back(facet_element);
            /// loop on every element connected to current facet and
            /// insert current facet in the first free spot of the
            /// subelement_to_element vector
            for (auto & loc_el : elements) {
              if (loc_el == ElementNull) {
                continue;
              }

              auto && subelements =
                  mesh_accessor.getSubelementToElement(loc_el);

              for (auto & el : subelements) {
                if (el != ElementNull) {
                  continue;
                }

                el = facet_element;
                break;
              }
            }

            /// reset connectivity in case a facet was found in
            /// between ghost and normal elements
            if (facet_ghost_type != ghost_type) {
              facet_ghost_type = ghost_type;
              connectivity_facets =
                  &mesh_accessor.getConnectivity(facet_type, facet_ghost_type);
              element_to_subelement = &mesh_accessor.getElementToSubelement(
                  facet_type, facet_ghost_type);
            }
          }
        }
      }
    }
  }

  mesh_facets.sendEvent(event);

  // restore the parent of mesh_facet
  if (mesh_facets_parent != nullptr) {
    mesh_facets.defineMeshParent(*mesh_facets_parent);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::renumberMeshNodes(Mesh & mesh,
                                  Array<UInt> & local_connectivities,
                                  UInt nb_local_element, UInt nb_ghost_element,
                                  ElementType type,
                                  Array<UInt> & old_nodes_numbers) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  std::map<UInt, UInt> renumbering_map;
  for (UInt i = 0; i < old_nodes_numbers.size(); ++i) {
    renumbering_map[old_nodes_numbers(i)] = i;
  }

  /// renumber the nodes
  renumberNodesInConnectivity(local_connectivities,
                              (nb_local_element + nb_ghost_element) *
                                  nb_nodes_per_element,
                              renumbering_map);

  old_nodes_numbers.resize(renumbering_map.size());
  for (auto & renumber_pair : renumbering_map) {
    old_nodes_numbers(renumber_pair.second) = renumber_pair.first;
  }
  renumbering_map.clear();

  MeshAccessor mesh_accessor(mesh);

  /// copy the renumbered connectivity to the right place
  auto & local_conn = mesh_accessor.getConnectivity(type);
  local_conn.resize(nb_local_element);

  if (nb_local_element > 0) {
    memcpy(local_conn.storage(), local_connectivities.storage(),
           nb_local_element * nb_nodes_per_element * sizeof(UInt));
  }

  auto & ghost_conn = mesh_accessor.getConnectivity(type, _ghost);
  ghost_conn.resize(nb_ghost_element);

  if (nb_ghost_element > 0) {
    std::memcpy(ghost_conn.storage(),
                local_connectivities.storage() +
                    nb_local_element * nb_nodes_per_element,
                nb_ghost_element * nb_nodes_per_element * sizeof(UInt));
  }

  auto & ghost_counter = mesh_accessor.getGhostsCounters(type, _ghost);
  ghost_counter.resize(nb_ghost_element, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::renumberNodesInConnectivity(
    Array<UInt> & list_nodes, UInt nb_nodes,
    std::map<UInt, UInt> & renumbering_map) {
  AKANTU_DEBUG_IN();

  UInt * connectivity = list_nodes.storage();
  UInt new_node_num = renumbering_map.size();
  for (UInt n = 0; n < nb_nodes; ++n, ++connectivity) {
    UInt & node = *connectivity;
    auto it = renumbering_map.find(node);
    if (it == renumbering_map.end()) {
      UInt old_node = node;
      renumbering_map[old_node] = new_node_num;
      node = new_node_num;
      ++new_node_num;
    } else {
      node = it->second;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::purifyMesh(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  std::map<UInt, UInt> renumbering_map;

  RemovedNodesEvent remove_nodes(mesh, AKANTU_CURRENT_FUNCTION);
  Array<UInt> & nodes_removed = remove_nodes.getList();

  for (auto ghost_type : ghost_types) {
    for (auto type :
         mesh.elementTypes(_all_dimensions, ghost_type, _ek_not_defined)) {
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
      UInt nb_element(connectivity.size());

      renumberNodesInConnectivity(
          connectivity, nb_element * nb_nodes_per_element, renumbering_map);
    }
  }

  Array<UInt> & new_numbering = remove_nodes.getNewNumbering();
  std::fill(new_numbering.begin(), new_numbering.end(), UInt(-1));

  for (auto && pair : renumbering_map) {
    new_numbering(std::get<0>(pair)) = std::get<1>(pair);
  }

  for (UInt i = 0; i < new_numbering.size(); ++i) {
    if (new_numbering(i) == UInt(-1)) {
      nodes_removed.push_back(i);
    }
  }

  mesh.sendEvent(remove_nodes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::flipFacets(
    Mesh & mesh_facets,
    const ElementTypeMapArray<UInt> & remote_global_connectivities,
    GhostType gt_facet) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();

  /// get global connectivity for local mesh
  ElementTypeMapArray<UInt> local_global_connectivities(
      "local_global_connectivity", mesh_facets.getID(),
      mesh_facets.getMemoryID());

  local_global_connectivities.initialize(
      mesh_facets, _spatial_dimension = spatial_dimension - 1,
      _ghost_type = gt_facet, _with_nb_nodes_per_element = true,
      _with_nb_element = true);

  mesh_facets.getGlobalConnectivity(local_global_connectivities);
  MeshAccessor mesh_accessor(mesh_facets);

  /// loop on every facet
  for (auto type_facet :
       mesh_facets.elementTypes(spatial_dimension - 1, gt_facet)) {

    auto & connectivity = mesh_accessor.getConnectivity(type_facet, gt_facet);
    auto & local_global_connectivity =
        local_global_connectivities(type_facet, gt_facet);
    const auto & remote_global_connectivity =
        remote_global_connectivities(type_facet, gt_facet);

    auto & element_per_facet =
        mesh_accessor.getElementToSubelementNC(type_facet, gt_facet);
    auto & subfacet_to_facet =
        mesh_accessor.getSubelementToElementNC(type_facet, gt_facet);

    auto nb_nodes_per_facet = connectivity.getNbComponent();
    auto nb_nodes_per_P1_facet =
        Mesh::getNbNodesPerElement(Mesh::getP1ElementType(type_facet));

    for (auto && data :
         zip(make_view(connectivity, nb_nodes_per_facet),
             make_view(local_global_connectivity, nb_nodes_per_facet),
             make_view(remote_global_connectivity, nb_nodes_per_facet),
             make_view(subfacet_to_facet, subfacet_to_facet.getNbComponent()),
             make_view(element_per_facet))) {

      auto & conn = std::get<0>(data);
      auto & local_gconn = std::get<1>(data);
      const auto & remote_gconn = std::get<2>(data);

      /// skip facet if connectivities are the same
      if (local_gconn == remote_gconn) {
        continue;
      }

      /// re-arrange connectivity
      auto conn_tmp = conn;
      auto begin = local_gconn.begin();
      auto end = local_gconn.end();

      AKANTU_DEBUG_ASSERT(std::is_permutation(begin, end, remote_gconn.begin()),
                          "This facets are not just permutation of each other, "
                              << local_gconn << " and " << remote_gconn);

      for (auto && data : enumerate(remote_gconn)) {
        auto it = std::find(begin, end, std::get<1>(data));
        AKANTU_DEBUG_ASSERT(it != end, "Node not found");
        UInt new_position = it - begin;
        conn(new_position) = conn_tmp(std::get<0>(data));
        ;
      }
      // std::transform(remote_gconn.begin(), remote_gconn.end(), conn.begin(),
      //                [&](auto && gnode) {
      //                  auto it = std::find(begin, end, gnode);
      //                  AKANTU_DEBUG_ASSERT(it != end, "Node not found");
      //                  return conn_tmp(it - begin);
      //                });

      /// if 3D, check if facets are just rotated
      if (spatial_dimension == 3) {
        auto begin = remote_gconn.begin();
        /// find first node
        auto it = std::find(begin, remote_gconn.end(), local_gconn(0));

        UInt n;
        UInt start = it - begin;
        /// count how many nodes in the received connectivity follow
        /// the same order of those in the local connectivity
        for (n = 1; n < nb_nodes_per_P1_facet &&
                    local_gconn(n) ==
                        remote_gconn((start + n) % nb_nodes_per_P1_facet);
             ++n) {
          ;
        }

        /// skip the facet inversion if facet is just rotated
        if (n == nb_nodes_per_P1_facet) {
          continue;
        }
      }

      /// update data to invert facet
      auto & element_per_facet = std::get<4>(data);
      if (element_per_facet[1] !=
          ElementNull) { // by convention the first facet
                         // cannot be a ElementNull
        std::swap(element_per_facet[0], element_per_facet[1]);
      }

      auto & subfacets_of_facet = std::get<3>(data);
      std::swap(subfacets_of_facet(0), subfacets_of_facet(1));
    }
  }

  AKANTU_DEBUG_OUT();
}

/*-------------------------------------------------------------------- */
Real MeshUtils::cosSharpAngleBetween2Facets(SolidMechanicsModel & model,
                                            const Element & facet1,
                                            const Element & facet2) {
  AKANTU_DEBUG_IN();

  auto & facets_fe_engine = model.getFEEngine("FacetsFEEngine");
  auto dim = model.getSpatialDimension();

  // normal to the first facet
  UInt nb_quad_points1 = facets_fe_engine.getNbIntegrationPoints(facet1.type);
  const auto & normals1 = facets_fe_engine.getNormalsOnIntegrationPoints(
      facet1.type, facet1.ghost_type);
  auto normal1_begin = normals1.begin(dim);
  Vector<Real> facet1_normal(normal1_begin[facet1.element * nb_quad_points1]);

  // normal to the second facet
  UInt nb_quad_points2 = facets_fe_engine.getNbIntegrationPoints(facet2.type);
  const auto & normals2 = facets_fe_engine.getNormalsOnIntegrationPoints(
      facet2.type, facet2.ghost_type);
  auto normal2_begin = normals2.begin(dim);
  Vector<Real> facet2_normal(normal2_begin[facet2.element * nb_quad_points2]);

  // get abs of dot product between two normals
  Real dot = facet1_normal.dot(facet2_normal);
  dot = std::abs(dot);

  AKANTU_DEBUG_OUT();
  return dot;
}

/*-------------------------------------------------------------------- */
Real MeshUtils::distanceBetween2Barycenters(const Mesh & mesh_facet1,
                                            const Mesh & mesh_facet2,
                                            const Element & facet1,
                                            const Element & facet2) {
  AKANTU_DEBUG_IN();
  auto dim = mesh_facet1.getSpatialDimension();
  dim = std::max(mesh_facet2.getSpatialDimension(), dim);

  Vector<Real> facet1_bary(dim);
  Vector<Real> facet2_bary(dim);
  mesh_facet1.getBarycenter(facet1, facet1_bary);
  mesh_facet2.getBarycenter(facet2, facet2_bary);
  auto dist = std::abs(facet1_bary.distance(facet2_bary));

  AKANTU_DEBUG_OUT();
  return dist;
}

/*-------------------------------------------------------------------- */
Real MeshUtils::distanceBetweenIncentersCorrected(const Mesh & mesh_facets,
                                                  const Element & facet1,
                                                  const Element & facet2) {
  AKANTU_DEBUG_IN();
  auto dim = mesh_facets.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim == 3,
                      "DistanceBetweenBarycentersCorrected works only in 3D");

  const Vector<Element> & subfacet_to_el1 =
      mesh_facets.getSubelementToElement(facet1);
  const Vector<Element> && subfacet_to_el2 =
      mesh_facets.getSubelementToElement(facet2);
  std::set<Element> intersection;

  for (auto & subfacet1 : subfacet_to_el1) {
    for (auto & subfacet2 : subfacet_to_el2) {
      if (subfacet1 == subfacet2) {
        intersection.emplace(subfacet1);
      }
    }
  }

  // verify if facets share 1 subfacet == connected but not coincide
  auto nb_shared_subfacets = intersection.size();
  AKANTU_DEBUG_ASSERT(
      nb_shared_subfacets == 1,
      "Facets " << facet1 << " and " << facet2
                << " are either not connected or coincide. Outputing zero");
  auto common_subfacet = *intersection.begin();
  auto & pos = mesh_facets.getNodes();
  const auto pos_it = pos.begin(dim);
  // auto subfacet_conn = mesh_facets.getConnectivity(common_subfacet.type,
  //                                                  common_subfacet.ghost_type);
  // auto nb_nodes_subfacet = subfacet_conn.getNbComponent();
  // const auto subfacet_nodes_it =
  //     make_view(subfacet_conn, nb_nodes_subfacet).begin();
  Vector<UInt> subfacet_nodes = mesh_facets.getConnectivity(common_subfacet);
  Vector<Real> A(dim);
  Vector<Real> AB(dim);
  A = Vector<Real>(pos_it[subfacet_nodes(0)]);
  AB = Vector<Real>(pos_it[subfacet_nodes(1)]) -
       Vector<Real>(pos_it[subfacet_nodes(0)]);

  Vector<Real> facet1_inc(dim);
  Vector<Real> facet2_inc(dim);
  mesh_facets.getIncenter(facet1, facet1_inc);
  mesh_facets.getIncenter(facet2, facet2_inc);

  Vector<Real> AInc1 = facet1_inc - A;
  Vector<Real> AInc2 = facet2_inc - A;
  Vector<Real> dif = (AInc1.norm() - AInc2.norm()) * AB / AB.norm();
  Vector<Real> facet2_inc_corrected = facet2_inc + dif;

  auto dist = std::abs(facet1_inc.distance(facet2_inc_corrected));

  AKANTU_DEBUG_OUT();
  return dist;
}
/*-------------------------------------------------------------------- */
Real MeshUtils::distanceBetweenIncenters(const Mesh & mesh_facets,
                                         const Element & facet1,
                                         const Element & facet2) {
  AKANTU_DEBUG_IN();
  auto dim = mesh_facets.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim == 3, "DistanceBetweenIncenters works only in 3D");

  Vector<Real> facet1_inc(dim);
  Vector<Real> facet2_inc(dim);
  mesh_facets.getIncenter(facet1, facet1_inc);
  mesh_facets.getIncenter(facet2, facet2_inc);

  auto dist = std::abs(facet1_inc.distance(facet2_inc));

  AKANTU_DEBUG_OUT();
  return dist;
}
/*-------------------------------------------------------------------- */
Real MeshUtils::getInscribedCircleDiameter(SolidMechanicsModel & model,
                                           const Element & el,
                                           bool initial_conf) {
  AKANTU_DEBUG_IN();
  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto & facets_fe_engine = model.getFEEngine("FacetsFEEngine");
  auto dim = mesh.getSpatialDimension();
  const auto & pos = mesh.getNodes();
  auto coordinates = pos;
  if (not initial_conf) {
    coordinates += model.getDisplacement();
  }

  auto nb_nodes_per_facet = mesh_facets.getNbNodesPerElement(el.type);
  Array<Real> coord(0, nb_nodes_per_facet * dim);
  Array<UInt> dummy_list(1, 1, el.element);
  facets_fe_engine.extractNodalToElementField(
      mesh_facets, coordinates, coord, el.type, el.ghost_type, dummy_list);
  Array<Real>::matrix_iterator coord_el = coord.begin(dim, nb_nodes_per_facet);
  Real facet_indiam = facets_fe_engine.getElementInradius(*coord_el, el.type);

  AKANTU_DEBUG_OUT();
  return facet_indiam;
}
/*-------------------------------------------------------------------- */
Real MeshUtils::getFacetArea(SolidMechanicsModel & model, const Element & el) {
  AKANTU_DEBUG_IN();

  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(el.type) == dim - 1,
                      "Element has a wrong dimension");
  AKANTU_DEBUG_ASSERT(Mesh::getKind(el.type) == _ek_regular,
                      "Element has a wrong kind: " << Mesh::getKind(el.type));

  auto & fe_engine = model.getFEEngine("FacetsFEEngine");
  Array<UInt> single_el_array;
  single_el_array.push_back(el.element);

  Array<Real> dummy(fe_engine.getNbIntegrationPoints(el.type), 1, 1.);
  Real el_area =
      fe_engine.integrate(dummy, el.type, el.ghost_type, single_el_array);

  AKANTU_DEBUG_OUT();
  return el_area;
}
/*-------------------------------------------------------------------- */
Real MeshUtils::getCurrentFacetArea(SolidMechanicsModel & model,
                                    const Element & el) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(el.type == _triangle_3,
                      "The only facet type supported is _triangle_3");
  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  const auto & pos = mesh.getNodes();
  auto coordinates = pos;
  coordinates += model.getDisplacement();
  auto coord_it = make_view(coordinates, dim).begin();

  auto && facet_conn = mesh_facets.getConnectivity(el.type, el.ghost_type);
  auto nb_nodes_facet = facet_conn.getNbComponent();
  auto facet_nodes_it = make_view(facet_conn, nb_nodes_facet).begin();
  auto facet_nodes = facet_nodes_it[el.element];

  // compute triangle's area from sides
  UInt A, B, C;
  Real a, b, c;
  A = facet_nodes(0);
  B = facet_nodes(1);
  C = facet_nodes(2);
  Vector<Real> AB = Vector<Real>(coord_it[B]) - Vector<Real>(coord_it[A]);
  Vector<Real> BC = Vector<Real>(coord_it[C]) - Vector<Real>(coord_it[B]);
  Vector<Real> AC = Vector<Real>(coord_it[C]) - Vector<Real>(coord_it[B]);
  a = AB.norm();
  b = BC.norm();
  c = AC.norm();
  auto s = (a + b + c) / 2;

  auto area = sqrt(s * (s - a) * (s - b) * (s - c));

  AKANTU_DEBUG_OUT();
  return area;
}
/*-------------------------------------------------------------------- */
std::pair<bool, Array<Element>>
MeshUtils::areFacetsConnected(const Mesh & mesh_facets, const Element & facet1,
                              const Element & facet2) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(
      facet1.type == facet2.type,
      "Different facet types are entered in areFacetsConnected function");

  Array<Element> shared_subs;
  bool facets_connected{false};

  const Vector<Element> & sub_to_facet1 =
      mesh_facets.getSubelementToElement(facet1);
  const Vector<Element> & sub_to_facet2 =
      mesh_facets.getSubelementToElement(facet2);
  Vector<Element> common_subfacets;

  for (UInt i : arange(sub_to_facet1.size())) {
    for (UInt j : arange(i, sub_to_facet2.size())) {
      auto sub1 = sub_to_facet1(i);
      auto sub2 = sub_to_facet2(j);
      if (sub1 == sub2) {
        shared_subs.push_back(sub1);
        facets_connected = true;
      }
    }
  }

  AKANTU_DEBUG_OUT();
  return std::make_pair(facets_connected, shared_subs);
}
/*-------------------------------------------------------------------- */
void MeshUtils::fillElementToSubElementsData(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  if (mesh.getNbElement(mesh.getSpatialDimension() - 1) == 0) {
    AKANTU_DEBUG_INFO("There are not facets, add them in the mesh file or call "
                      "the buildFacet method.");
    return;
  }

  UInt spatial_dimension = mesh.getSpatialDimension();
  ElementTypeMapArray<Real> barycenters("barycenter_tmp", mesh.getID(),
                                        mesh.getMemoryID());
  barycenters.initialize(mesh, _nb_component = spatial_dimension,
                         _spatial_dimension = _all_dimensions);

  Element element;
  for (auto ghost_type : ghost_types) {
    element.ghost_type = ghost_type;
    for (const auto & type : mesh.elementTypes(_all_dimensions, ghost_type)) {
      element.type = type;

      UInt nb_element = mesh.getNbElement(type, ghost_type);
      Array<Real> & barycenters_arr = barycenters(type, ghost_type);
      barycenters_arr.resize(nb_element);
      auto bary = barycenters_arr.begin(spatial_dimension);
      auto bary_end = barycenters_arr.end(spatial_dimension);

      for (UInt el = 0; bary != bary_end; ++bary, ++el) {
        element.element = el;
        mesh.getBarycenter(element, *bary);
      }
    }
  }

  MeshAccessor mesh_accessor(mesh);
  for (Int sp(spatial_dimension); sp >= 1; --sp) {
    if (mesh.getNbElement(sp) == 0) {
      continue;
    }

    for (auto ghost_type : ghost_types) {
      for (auto & type : mesh.elementTypes(sp, ghost_type)) {
        auto & subelement_to_element =
            mesh_accessor.getSubelementToElement(type, ghost_type);
        subelement_to_element.resize(mesh.getNbElement(type, ghost_type));
        subelement_to_element.set(ElementNull);
      }

      for (auto & type : mesh.elementTypes(sp - 1, ghost_type)) {
        auto & element_to_subelement =
            mesh_accessor.getElementToSubelement(type, ghost_type);
        element_to_subelement.resize(mesh.getNbElement(type, ghost_type));
        element_to_subelement.clear();
      }
    }

    CSR<Element> nodes_to_elements;
    buildNode2Elements(mesh, nodes_to_elements, sp);

    Element facet_element;

    for (auto ghost_type : ghost_types) {
      facet_element.ghost_type = ghost_type;
      for (const auto & type : mesh.elementTypes(sp - 1, ghost_type)) {
        facet_element.type = type;

        auto & element_to_subelement =
            mesh_accessor.getElementToSubelement(type, ghost_type);

        const auto & connectivity = mesh.getConnectivity(type, ghost_type);
        // element_to_subelement.resize(connectivity.size());

        for (auto && data : enumerate(
                 make_view(connectivity, mesh.getNbNodesPerElement(type)))) {
          const auto & facet = std::get<1>(data);
          facet_element.element = std::get<0>(data);

          std::map<Element, UInt> element_seen_counter;
          auto nb_nodes_per_facet =
              mesh.getNbNodesPerElement(Mesh::getP1ElementType(type));

          // count the number of node in common between the facet and the
          // other element connected to the nodes of the facet
          for (auto node : arange(nb_nodes_per_facet)) {
            for (auto & elem : nodes_to_elements.getRow(facet(node))) {
              auto cit = element_seen_counter.find(elem);
              if (cit != element_seen_counter.end()) {
                cit->second++;
              } else {
                element_seen_counter[elem] = 1;
              }
            }
          }

          // check which are the connected elements
          std::vector<Element> connected_elements;
          for (auto && cit : element_seen_counter) {
            if (cit.second == nb_nodes_per_facet) {
              connected_elements.push_back(cit.first);
            }
          }

          // add the connected elements as sub-elements
          for (auto & connected_element : connected_elements) {
            element_to_subelement(facet_element.element)
                .push_back(connected_element);
          }

          // add the element as sub-element to the connected elements
          for (auto & connected_element : connected_elements) {
            Vector<Element> subelements_to_element =
                mesh.getSubelementToElement(connected_element);

            // find the position where to insert the element
            auto it = std::find(subelements_to_element.begin(),
                                subelements_to_element.end(), ElementNull);

            AKANTU_DEBUG_ASSERT(
                it != subelements_to_element.end(),
                "The element "
                    << connected_element << " seems to have too many facets!! ("
                    << (it - subelements_to_element.begin()) << " < "
                    << mesh.getNbFacetsPerElement(connected_element.type)
                    << ")");

            *it = facet_element;
          }
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
