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

#ifndef AKANTU_MESH_SPHERE_INTERSECTOR_HH_
#define AKANTU_MESH_SPHERE_INTERSECTOR_HH_

#include "aka_common.hh"
#include "mesh_geom_intersector.hh"

#include "mesh_geom_common.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

template <Int dim, ElementType type>
class MeshSphereIntersector
    : public MeshGeomIntersector<dim, type, Line_arc<cgal::Spherical>,
                                 cgal::Spherical::Sphere_3, cgal::Spherical> {
  using SK = cgal::Spherical;
  using K = cgal::Cartesian;

  /// Parent class type
  using parent_type =
      MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK>;

  /// Result of intersection function type
  using result_type =
      typename IntersectionTypeHelper<TreeTypeHelper<Triangle<K>, K>,
                                      K::Segment_3>::intersection_type;

  /// Pair of intersection points and element id
  using pair_type = std::pair<SK::Circular_arc_point_3, Int>;

public:
  /// Construct from mesh
  explicit MeshSphereIntersector(Mesh & mesh);

  /// Destructor
  virtual ~MeshSphereIntersector();

public:
  /// Construct the primitive tree object
  void constructData(GhostType ghost_type = _not_ghost) override;

  /**
   * @brief Computes the intersection of the mesh with a sphere
   */
  void computeIntersectionQuery(const SK::Sphere_3 & /* query */) override {
    AKANTU_ERROR("This function is not implemented for spheres (It was "
                 "to generic and has been replaced by "
                 "computeMeshQueryIntersectionPoint");
  }

  /**
   * Compute intersection points between the mesh primitives (segments) and a
   * query (surface in 3D or a curve in 2D), double intersection points for the
   * same primitives are not considered. A maximum is set to the number of
   * intersection nodes per element: 2 in 2D and 4 in 3D
   */
  void computeMeshQueryIntersectionPoint(const SK::Sphere_3 & query,
                                         Int nb_old_nodes) override;

  /// Build the IGFEM mesh
  void
  buildResultFromQueryList(const std::list<SK::Sphere_3> & /*query*/) override {
    AKANTU_ERROR("This function is no longer implemented to split "
                 "geometrical operations and dedicated result "
                 "construction");
  }

  /// Set the tolerance
  void setToleranceIntersectionOnNode(UInt tol) {
    this->tol_intersection_on_node = tol;
  }

protected:
  /// tolerance for which the intersection is considered on the mesh node
  /// (relative to the segment lenght)
  Real tol_intersection_on_node;
};

} // namespace akantu

#include "mesh_sphere_intersector_tmpl.hh"

#endif // AKANTU_MESH_SPHERE_INTERSECTOR_HH_
