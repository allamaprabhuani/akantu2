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
#include "aka_common.hh"
#include "aka_csr.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MESH_UTILS_HH_
#define AKANTU_MESH_UTILS_HH_

namespace akantu {

class MeshUtils {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// build a CSR<Element> that contains for each node the list of connected
  /// elements of a given spatial dimension
  AKANTU_EXPORT static void
  buildNode2Elements(const Mesh & mesh, CSR<Element> & node_to_elem,
                     Int spatial_dimension = _all_dimensions);

  /// build a CSR<UInt> that contains for each node the number of
  /// the connected elements of a given ElementType
  static void
  buildNode2ElementsElementTypeMap(const Mesh & mesh, CSR<Idx> & node_to_elem,
                                   ElementType type,
                                   GhostType ghost_type = _not_ghost);

  /// build the facets elements on the boundaries of a mesh
  AKANTU_EXPORT static void buildFacets(Mesh & mesh);

  /// build all the facets elements: boundary and internals and store them in
  /// the mesh_facets for element of dimension from_dimension to to_dimension
  AKANTU_EXPORT static void buildAllFacets(const Mesh & mesh,
                                           Mesh & mesh_facets,
                                           Int from_dimension,
                                           Int to_dimension);

  /// build all the facets elements: boundary and internals and store them in
  /// the mesh_facets
  AKANTU_EXPORT static void
  buildAllFacets(const Mesh & mesh, Mesh & mesh_facets, Int to_dimension = 0);

  /// build facets for a given spatial dimension
  static void buildFacetsDimension(const Mesh & mesh, Mesh & mesh_facets,
                                   bool boundary_only, Int dimension);

  /// take the local_connectivity array as the array of local and ghost
  /// connectivity, renumber the nodes and set the connectivity of the mesh
  static void renumberMeshNodes(Mesh & mesh, Array<Idx> & local_connectivities,
                                Int nb_local_element, Int nb_ghost_element,
                                ElementType type, Array<Idx> & old_nodes);

  /// compute pbc pair for a given direction
  static void computePBCMap(const Mesh & mymesh, Int dir,
                            std::map<Idx, Idx> & pbc_pair);
  /// compute pbc pair for a surface pair
  static void computePBCMap(const Mesh & mymesh,
                            const std::pair<ID, ID> & surface_pair,
                            std::map<Idx, Idx> & pbc_pair);

  /// remove not connected nodes /!\ this functions renumbers the nodes.
  AKANTU_EXPORT static void purifyMesh(Mesh & mesh);

  /// fill the subelement to element and the elements to subelements data
  static void fillElementToSubElementsData(Mesh & mesh);

  /// flip facets based on global connectivity
  static void
  flipFacets(Mesh & mesh_facets,
             const ElementTypeMapArray<Idx> & remote_global_connectivities,
             GhostType gt_facet);

private:
  /// match pairs that are on the associated pbc's
  static void matchPBCPairs(const Mesh & mymesh, Int dir,
                            Array<Idx> & selected_left,
                            Array<Idx> & selected_right,
                            std::map<Idx, Idx> & pbc_pair);

  /// function used by all the renumbering functions
  static void renumberNodesInConnectivity(Array<Idx> & list_nodes, Int nb_nodes,
                                          std::map<Idx, Idx> & renumbering_map);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "mesh_utils_inline_impl.hh"

#endif /* AKANTU_MESH_UTILS_HH_ */
