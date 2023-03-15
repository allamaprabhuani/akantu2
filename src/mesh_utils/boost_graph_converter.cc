#include "boost_graph_converter.hh"
#include "mesh.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <boost/graph/graphviz.hpp>
#include <map>
/* -------------------------------------------------------------------------- */

namespace akantu {

BoostGraphConverter::BoostGraphConverter(const Mesh & mesh) {
  // adding all vertices
  auto add_vertex = [&](auto element) {
    auto && vertex = boost::add_vertex(this->graph);
    element_to_vertex[element] = vertex;
    vertex_to_element[vertex] = element;
  };

  //  add_vertex(ElementNull);

  auto && mesh_facets = mesh.getMeshFacets();

  for_each_element(
      mesh, [&](auto && element) { add_vertex(element); },
      _element_kind = _ek_regular);
  for_each_element(
      mesh_facets, [&](auto && element) { add_vertex(element); },
      _spatial_dimension = _all_dimensions);
  for_each_element(
      mesh, [&](auto && element) { add_vertex(element); },
      _element_kind = _ek_cohesive);

  auto get_vertex = [&](auto && element) {
    if (element == ElementNull) {
      return -1;
    }

    auto it = element_to_vertex.find(element);
    if (it == element_to_vertex.end()) {
      AKANTU_EXCEPTION("Element " << element << " not in element to vertex");
    }
    return it->second;
  };

  auto add_edge = [&](auto && element1, auto && element2,
                      EdgeType type = _none) {
    auto vertex1 = get_vertex(element1);
    auto vertex2 = get_vertex(element2);

    if (vertex1 == -1 or vertex2 == -1) {
      return;
    }

    // std::cout << element1 << " -> " << element2 << " " << type << std::endl;
    auto && data = boost::add_edge(vertex1, vertex2, graph);
    edge_type[data.first] = type;
  };

  auto add_neighbors = [&](auto && element) {
    const auto & element_facets =
        mesh_facets.getSubelementToElement().get(element);
    for (auto && facet : element_facets) {
      add_edge(element, facet, _up);
      const auto & connected_elements_to_facets =
          const_cast<const Mesh &>(mesh_facets).getElementToSubelement(facet);
      for (auto && celement : connected_elements_to_facets) {
        add_edge(facet, celement, _down);
      }
    }
  };
  // adding edges
  for_each_element(mesh, add_neighbors, _element_kind = _ek_regular);
  for (auto s : arange(1, mesh.getSpatialDimension())) {
    for_each_element(mesh_facets, add_neighbors, _spatial_dimension = s);
  }
  for_each_element(mesh, add_neighbors, _element_kind = _ek_cohesive);
}

void BoostGraphConverter::toGraphviz(const std::string & filename) const {
  std::ofstream fout;
  fout.open(filename);

  boost::write_graphviz(
      fout, graph,
      [&](std::ostream & out, vertex_descriptor vertex) {
        auto it = vertex_to_element.find(vertex);
        if (it == vertex_to_element.end()) {
          return;
        }

        auto && element = it->second;
        std::array<std::string, 4> colors{"gold", "cadetblue1", "coral1",
                                          "antiquewhite2"};

        std::string color{"deeppink4"}, bgcolor{"deeppink4"};
        if (element != ElementNull) {
          color = bgcolor = colors[Mesh::getSpatialDimension(element.type)];
        }

        out << "[label=\"" << std::to_string(element)
            << "\", bgcolor=" << bgcolor << ", color=" << color << "]";
      },
      [&](std::ostream & out, edge_descriptor edge) {
        auto && type = edge_type.at(edge);
        std::map<EdgeType, std::string> colors{
            {_none, "black"}, {_up, "darkorange"}, {_down, "blue2"}};
        out << "[color=" << colors[type] << "]";
      });
}

} // namespace akantu
