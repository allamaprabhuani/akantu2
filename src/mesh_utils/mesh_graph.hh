/**
 * @file   mesh_graph.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Fri Jun 06 2014
 *
 * @brief  graph for the identification of mesh fragments
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_GRAPH_HH__
#define __AKANTU_MESH_GRAPH_HH__

#include <iostream>

#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
#elif (defined(__GNUC__) || defined(__GNUG__))
#  define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#  if GCC_VERSION > 40800
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#  endif
#endif


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "mesh.hh"
#include "element_type_map.hh"

#define DEBUG_GRAPH 1

__BEGIN_AKANTU__

using std::cout;
using std::endl;


class MeshGraph {

public:
  
  // graph type definitions
  typedef boost::property<boost::vertex_index1_t, UInt> vertex_property;
  typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, vertex_property> graph_type;
  typedef boost::graph_traits<graph_type> graph_traits;
  typedef graph_traits::vertex_iterator vertex_iter;
  typedef graph_traits::vertex_descriptor vertex_type;


  MeshGraph(Mesh&, Mesh&);

//  template <class filter_type>
//  MeshGraph(Mesh&, Mesh&, filter_type& = nullptr);

  const Array<UInt>& components(ElementType type, GhostType ghost_type = _not_ghost) const
  { return components_(type,ghost_type); }
  
  size_t numComponents() const
  { return num_; }
  
  
  friend std::ostream& operator<<(std::ostream& os, const MeshGraph& mg) {
    
    os<<"Vertices: ";
    std::pair<vertex_iter, vertex_iter> vp;
    for (vp = boost::vertices(mg.g_); vp.first != vp.second; ++vp.first)
      os << *vp.first << " ";
    os << std::endl;
    
    // Iterate through the edges and print them out
    os<<"Edges: ";
    typedef graph_traits::edge_iterator edge_iter;
    std::pair<edge_iter, edge_iter> ep;
    edge_iter ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(mg.g_); ei != ei_end; ++ei)
      os << *ei<<" ";
    os<<endl;
    
    return os;
  }
  
private:

  void connectedComponents(Mesh& mesh) {
    
    // initialize container
    mesh.initElementTypeMapArray(components_,1,mesh.getSpatialDimension(),false,_ek_regular,true);

    // obtain a property map for the vertex_index property
    boost::property_map<graph_type, boost::vertex_index1_t>::type index = boost::get(boost::vertex_index1_t(), g_);

    // call boost::connected_components
    std::vector<int> c(boost::num_vertices(g_));
    int num = boost::connected_components(g_, &c[0]);
    
#ifdef DEBUG_GRAPH
    
    std::vector<int>::size_type i;
    cout << "Total number of components: " << num << endl;
    for (i = 0; i != c.size(); ++i)
      cout << "Vertex " << i <<" is in component " << c[i] << endl;
    cout << endl;
#endif
    
    // fill ElementTypeMap components array
    std::pair<vertex_iter, vertex_iter> vp;
    size_t k=0;
    std::set<size_t> comp;
    for (vp = boost::vertices(g_); vp.first != vp.second; ++vp.first) {
     
      vertex_type u = *vp.first;
      UInt id = index[u];
      Element el = mesh.linearizedToElement(id);
      Array<UInt> &array = components_(el.type);
      array(el.element) = c[k];
      comp.insert(c[k++]);
    }
    num_ = comp.size();
  }

  
  graph_type g_;                       //!< Mesh graph
  ElementTypeMapArray<UInt> components_;       //!< Array that stores the components
  size_t num_;                         //!< Number of components
};


#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
#elif defined(__GNUG__)
#  if GCC_VERSION > 40800
#    pragma GCC diagnostic pop
#  endif
#endif


__END_AKANTU__

#endif /* __AKANTU_MESH_GRAPH_HH__ */
