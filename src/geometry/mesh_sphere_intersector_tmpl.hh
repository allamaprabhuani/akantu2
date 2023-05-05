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

#ifndef AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH_
#define AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH_

#include "aka_common.hh"
#include "mesh_geom_common.hh"
#include "mesh_sphere_intersector.hh"
#include "tree_type_helper.hh"

namespace akantu {

template <Int dim, ElementType type>
MeshSphereIntersector<dim, type>::MeshSphereIntersector(Mesh & mesh)
    : parent_type(mesh), tol_intersection_on_node(1e-10) {
#if defined(AKANTU_IGFEM)
  if constexpr ((type == _triangle_3) or (type == _igfem_triangle_4) or
                (type == _igfem_triangle_5)) {
    this->nb_seg_by_el = 3;
  } else {
    AKANTU_ERROR("Not ready for mesh type " << type);
  }
#else
  if (type != _triangle_3) {
    AKANTU_ERROR("Not ready for mesh type " << type);
  }
#endif

  // initialize the intersection pointsss array with the spatial dimension
  this->intersection_points = new Array<Real>(0, dim);
  //  A maximum is set to the number of intersection nodes per element to limit
  //  the size of new_node_per_elem: 2 in 2D and 4 in 3D
  this->new_node_per_elem = new Array<UInt>(0, 1 + 4 * (dim - 1));
}

template <Int dim, ElementType type>
MeshSphereIntersector<dim, type>::~MeshSphereIntersector() {
  delete this->new_node_per_elem;
  delete this->intersection_points;
}

template <Int dim, ElementType type>
void MeshSphereIntersector<dim, type>::constructData(GhostType ghost_type) {

  this->new_node_per_elem->resize(this->mesh.getNbElement(type, ghost_type));
  this->new_node_per_elem->clear();

  MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK>::constructData(
      ghost_type);
}

template <Int dim, ElementType type>
void MeshSphereIntersector<dim, type>::computeMeshQueryIntersectionPoint(
    const SK::Sphere_3 & query, Int nb_old_nodes) {
  /// function to replace computeIntersectionQuery in a more generic geometry
  /// module version
  // The newNodeEvent is not send from this method who only compute the
  // intersection points
  AKANTU_DEBUG_IN();

  auto & nodes = this->mesh.getNodes();
  auto nb_node = nodes.size() + this->intersection_points->size();

  // Tolerance for proximity checks should be defined by user
  Real global_tolerance = Math::getTolerance();
  Math::setTolerance(tol_intersection_on_node);
  using sk_inter_res = boost::variant<pair_type>;

  // loop on the primitives (segments)
  for (auto && primitive : this->factory.getPrimitiveList()) {
    std::list<sk_inter_res> s_results;
    CGAL::intersection(primitive, query, std::back_inserter(s_results));

    if (s_results.size() != 1) { // not just one point
      continue;
    }

    auto * pair = boost::get<pair_type>(&s_results.front());

    if (not pair and pair->second != 1) {
      continue;
    }

    // the intersection point written as a vector
    Vector<Real> new_node(dim, 0.0);
    cgal::Cartesian::Point_3 point(CGAL::to_double(pair->first.x()),
                                   CGAL::to_double(pair->first.y()),
                                   CGAL::to_double(pair->first.z()));
    for (Int i = 0; i < dim; i++) {
      new_node(i) = point[i];
    }

    /// boolean to decide wheter intersection point is on a standard node
    /// of the mesh or not
    bool is_on_mesh = false;
    /// boolean to decide if this intersection point has been already
    /// computed for a neighbor element
    bool is_new = true;

    // check if intersection point has already been computed
    auto old_node_id = nb_old_nodes;
    // check if we already compute this intersection and add it as a node
    // for a neighbor element of another type
    auto node_view = make_view<dim>(nodes);
    // loop on the nodes from nb_old_nodes
    for (auto && existing_node :
         range(node_view.begin() + nb_old_nodes, node_view.end())) {
      if (new_node.isApprox(existing_node)) {
        is_new = false;
        break;
      }
      ++old_node_id;
    }

    if (is_new) {
      for (auto && intersection_point :
           make_view<dim>(*this->intersection_points)) {
        if (new_node.isApprox(intersection_point)) {
          is_new = false;
          break;
        }
      }
    }

    // get the initial and final points of the primitive (segment) and
    // write them as vectors
    cgal::Cartesian::Point_3 source_cgal(
        CGAL::to_double(primitive.source().x()),
        CGAL::to_double(primitive.source().y()),
        CGAL::to_double(primitive.source().z()));
    cgal::Cartesian::Point_3 target_cgal(
        CGAL::to_double(primitive.target().x()),
        CGAL::to_double(primitive.target().y()),
        CGAL::to_double(primitive.target().z()));

    Vector<Real, dim> source, target;
    for (Int i = 0; i < dim; i++) {
      source(i) = source_cgal[i];
      target(i) = target_cgal[i];
    }

    // Check if we are close from a node of the primitive (segment)
    if (source.isApprox(new_node) or target.isApprox(new_node)) {
      is_on_mesh = true;
      is_new = false;
    }

    // if the intersection point is a new one add it to the
    // list
    if (is_new) {
      this->intersection_points->push_back(new_node);
      nb_node++;
    }

    // deduce the element id
    auto element_id = primitive.id();

    // fill the new_node_per_elem array
    if (not is_on_mesh) { // if the node is not on a mesh node
      auto & nb_new_nodes_per_el = (*this->new_node_per_elem)(element_id, 0);
      nb_new_nodes_per_el += 1;
      AKANTU_DEBUG_ASSERT(
          2 * nb_new_nodes_per_el < this->new_node_per_elem->getNbComponent(),
          "You might have to interface crossing the same material");
      (*this->new_node_per_elem)(element_id, (2 * nb_new_nodes_per_el) - 1) =
          old_node_id;
      (*this->new_node_per_elem)(element_id, 2 * nb_new_nodes_per_el) =
          primitive.segId();
    }
  }

  Math::setTolerance(global_tolerance);

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif // AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH_
