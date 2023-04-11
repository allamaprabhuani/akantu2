/**
 * Copyright (©) 2022-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_BOOST_GRAPH_CONVERTER_HH_
#define AKANTU_BOOST_GRAPH_CONVERTER_HH_

/* -------------------------------------------------------------------------- */
#include "element.hh"
/* -------------------------------------------------------------------------- */
#include <boost/graph/adjacency_list.hpp>
#include <map>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {
class Mesh;
} // namespace akantu

namespace akantu {

class BoostGraphConverter {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  BoostGraphConverter(const Mesh & mesh);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void toGraphviz(const std::string & filename) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  enum EdgeType {
    _none,
    _up,
    _down,
  };

  using Graph = boost::adjacency_list<>;

  using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
  using edge_descriptor = boost::graph_traits<Graph>::edge_descriptor;

  using VertexInfo = std::map<vertex_descriptor, Element>;
  using EdgeInfo = std::map<edge_descriptor, EdgeType>;

  Graph graph;
  std::map<Element, Int> element_to_vertex;

  VertexInfo vertex_to_element;
  EdgeInfo edge_type;
};

} // namespace akantu

#endif // AKANTU_BOOST_GRAPH_CONVERTER_HH_
