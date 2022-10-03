/**
 * @file   embedded_interface_intersector.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri May 01 2015
 * @date last modification: Tue May 21 2019
 *
 * @brief  Class that loads the interface from mesh and computes intersections
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "embedded_interface_intersector.hh"
#include "mesh_segment_intersector.hh"

/// Helper macro for types in the mesh. Creates an intersector and computes
/// intersection queries
#define INTERFACE_INTERSECTOR_CASE(dim, type)                                  \
  do {                                                                         \
    MeshSegmentIntersector<dim, type> intersector(this->mesh, interface_mesh); \
    for (auto && [name, segment_list] : name_to_primitives) {                  \
      intersector.setPhysicalName(name);                                       \
      intersector.buildResultFromQueryList(segment_list);                      \
    }                                                                          \
  } while (0)

#define INTERFACE_INTERSECTOR_CASE_2D(type) INTERFACE_INTERSECTOR_CASE(2, type)
#define INTERFACE_INTERSECTOR_CASE_3D(type) INTERFACE_INTERSECTOR_CASE(3, type)

namespace akantu {

EmbeddedInterfaceIntersector::EmbeddedInterfaceIntersector(
    Mesh & mesh, const Mesh & primitive_mesh)
    : MeshGeomAbstract(mesh),
      interface_mesh(mesh.getSpatialDimension(), "interface_mesh"),
      primitive_mesh(primitive_mesh) {
  // Initiating mesh connectivity and data
  interface_mesh.addConnectivityType(_segment_2, _not_ghost);
  interface_mesh.addConnectivityType(_segment_2, _ghost);
  interface_mesh.getElementalData<Element>("associated_element")
      .alloc(0, 1, _segment_2);
  interface_mesh.getElementalData<std::string>("physical_names")
      .alloc(0, 1, _segment_2);
}

void EmbeddedInterfaceIntersector::constructData(GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  const Int dim = this->mesh.getSpatialDimension();

  if (dim == 1) {
    AKANTU_ERROR(
        "No embedded model in 1D. Deactivate intersection initialization");
  }

  Array<std::string> * physical_names = nullptr;

  try {
    physical_names = &const_cast<Array<std::string> &>(
        this->primitive_mesh.getData<std::string>("physical_names",
                                                  _segment_2));
  } catch (debug::Exception & e) {
    AKANTU_ERROR("You must define physical names to reinforcements in "
                 "order to use the embedded model");
    throw e;
  }

  constexpr Int nb_nodes_per_element = Mesh::getNbNodesPerElement(_segment_2);
  auto connectivity =
      primitive_mesh.getConnectivity(_segment_2).begin(nb_nodes_per_element);

  std::map<std::string, std::list<K::Segment_3>> name_to_primitives;

  // Loop over the physical names and register segment lists in
  // name_to_primitives_map
  for (auto && [element_id, name] : enumerate(*physical_names)) {
    auto && el_connectivity = connectivity[element_id];

    auto segment = this->createSegment(el_connectivity);
    name_to_primitives[name].push_back(segment);
  }

  // Loop over the background types of the mesh
  for (auto type : this->mesh.elementTypes(dim, _not_ghost)) {
    // Used in AKANTU_BOOST_ELEMENT_SWITCH
    AKANTU_DEBUG_INFO("Computing intersections with background element type "
                      << type);

    tuple_dispatch<
        std::tuple<_element_type_triangle_3, _element_type_triangle_6,
                   _element_type_tetrahedron_4>>(
        [&](auto && enum_type) {
          constexpr auto type = std::decay_t<decltype(enum_type)>::value;
          constexpr auto dim = Mesh::getSpatialDimension(type);
          MeshSegmentIntersector<dim, type> intersector(this->mesh,
                                                        interface_mesh);
          for (auto && [name, segment_list] : name_to_primitives) {
            intersector.setPhysicalName(name);
            intersector.buildResultFromQueryList(segment_list);
          }
        },
        type);
  }

  AKANTU_DEBUG_OUT();
}

K::Segment_3
EmbeddedInterfaceIntersector::createSegment(const Vector<Idx> & connectivity) {
  AKANTU_DEBUG_IN();

  std::unique_ptr<K::Point_3> source;
  std::unique_ptr<K::Point_3> target;
  const auto & nodes = this->primitive_mesh.getNodes();

  if (this->mesh.getSpatialDimension() == 2) {
    source = std::make_unique<K::Point_3>(nodes(connectivity(0), 0),
                                          nodes(connectivity(0), 1), 0.);
    target = std::make_unique<K::Point_3>(nodes(connectivity(1), 0),
                                          nodes(connectivity(1), 1), 0.);
  } else if (this->mesh.getSpatialDimension() == 3) {
    source = std::make_unique<K::Point_3>(nodes(connectivity(0), 0),
                                          nodes(connectivity(0), 1),
                                          nodes(connectivity(0), 2));
    target = std::make_unique<K::Point_3>(nodes(connectivity(1), 0),
                                          nodes(connectivity(1), 1),
                                          nodes(connectivity(1), 2));
  }

  K::Segment_3 segment(*source, *target);

  AKANTU_DEBUG_OUT();
  return segment;
}

} // namespace akantu

#undef INTERFACE_INTERSECTOR_CASE
#undef INTERFACE_INTERSECTOR_CASE_2D
#undef INTERFACE_INTERSECTOR_CASE_3D
