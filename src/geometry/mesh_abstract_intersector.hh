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

#ifndef AKANTU_MESH_ABSTRACT_INTERSECTOR_HH_
#define AKANTU_MESH_ABSTRACT_INTERSECTOR_HH_

#include "aka_common.hh"
#include "mesh_geom_abstract.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * @brief Class used to perform intersections on a mesh and construct output
 * data
 */
template <class Query> class MeshAbstractIntersector : public MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshAbstractIntersector(Mesh & mesh);

  /// Destructor
  ~MeshAbstractIntersector() override = default;

public:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the new_node_per_elem array
  AKANTU_GET_MACRO(NewNodePerElem, *new_node_per_elem, const Array<Int> &);
  /// get the intersection_points array
  AKANTU_GET_MACRO(IntersectionPoints, *intersection_points,
                   const Array<Real> &);
  /// get the nb_seg_by_el UInt
  AKANTU_GET_MACRO(NbSegByEl, nb_seg_by_el, Int);

  /**
   * @brief Compute the intersection with a query object
   *
   * This function needs to be implemented for every subclass. It computes the
   * intersections
   * with the tree of primitives and creates the data for the user.
   *
   * @param query the CGAL primitive of the query object
   */
  virtual void computeIntersectionQuery(const Query & query) = 0;

  /// Compute intersection points between the mesh primitives (segments) and a
  /// query (surface in 3D or a curve in 2D), double intersection points for the
  /// same primitives are not considered. A maximum intersection node per
  /// element is set : 2 in 2D and 4 in 3D
  virtual void computeMeshQueryIntersectionPoint(const Query & query,
                                                 Int nb_old_nodes) = 0;

  /// Compute intersection between the mesh and a list of queries
  virtual void
  computeIntersectionQueryList(const std::list<Query> & query_list);

  /// Compute intersection points between the mesh and a list of queries
  virtual void
  computeMeshQueryListIntersectionPoint(const std::list<Query> & query_list,
                                        Int nb_old_nodes);

  /// Compute whatever result is needed from the user (should be move to the
  /// appropriate specific classe for genericity)
  virtual void
  buildResultFromQueryList(const std::list<Query> & query_list) = 0;

protected:
  /// new node per element (column 0: number of new nodes, then odd is the
  /// intersection node number and even the ID of the intersected segment)
  std::unique_ptr<Array<Int>> new_node_per_elem;

  /// intersection output: new intersection points
  /// (computeMeshQueryListIntersectionPoint)
  std::unique_ptr<Array<Real>> intersection_points;

  /// number of segment in a considered element of the templated type of element
  /// specialized intersector
  const Int nb_seg_by_el{0};
};

} // namespace akantu

#include "mesh_abstract_intersector_tmpl.hh"

#endif // AKANTU_MESH_ABSTRACT_INTERSECTOR_HH_
