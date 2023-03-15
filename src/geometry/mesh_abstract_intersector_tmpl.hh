/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH_
#define AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH_

#include "aka_common.hh"
#include "mesh_abstract_intersector.hh"

namespace akantu {

template <class Query>
MeshAbstractIntersector<Query>::MeshAbstractIntersector(Mesh & mesh)
    : MeshGeomAbstract(mesh) {}

template <class Query>
void MeshAbstractIntersector<Query>::computeIntersectionQueryList(
    const std::list<Query> & query_list) {
  for (auto && query : query_list) {
    computeIntersectionQuery(query);
  }
}

template <class Query>
void MeshAbstractIntersector<Query>::computeMeshQueryListIntersectionPoint(
    const std::list<Query> & query_list, Int nb_old_nodes) {
  for (auto && query : query_list) {
    computeMeshQueryIntersectionPoint(query, nb_old_nodes);
  }
}

} // namespace akantu

#endif // AKANTU_MESH_ABSTRACT_INTERSECTOR_TMPL_HH_
