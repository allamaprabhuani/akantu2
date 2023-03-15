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
