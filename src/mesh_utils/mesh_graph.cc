/**
 * @file   mesh_graph.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Thu Jun 05 2014
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

#include "mesh_graph.hh"
#include "mesh_utils.hh"


using namespace boost;

__BEGIN_AKANTU__


MeshGraph::MeshGraph(Mesh& mesh, Mesh& facets) : g_(), components_("components","mesh_graph") {
  
  typedef std::map<UInt, vertex_type> check_map;
  typedef typename check_map::iterator check_iterator;
  
  // get mesh
  UInt d = mesh.getSpatialDimension();
  
  // obtain a property map for the vertex_index property
  boost::property_map<graph_type, boost::vertex_index1_t>::type index = boost::get(boost::vertex_index1_t(), g_);
  
  // container to keep track of inserted vertices
  check_map inserted;
  
  // loop over facets
  for (Mesh::type_iterator it_f = facets.firstType(d - 1); it_f != facets.lastType(d - 1); ++it_f) {
    
    UInt nb_facet = facets.getNbElement(*it_f);
    
    /// get the element connected to a subelement
    const Array< std::vector<Element> > &facet_elems = facets.getElementToSubelement(*it_f);
    
    // loop over facets to construct mesh dual graph
    for (UInt i=0; i<nb_facet; ++i) {
      
      // skip if facet is on the boundary
      if (facet_elems(i)[0].ghost_type == _ghost || facet_elems(i)[1].ghost_type == _ghost)
        continue;
      
      std::vector<vertex_type> u;
      for (size_t j=0; j<2; ++j) {
        
        if (facet_elems(i)[j] == ElementNull)
          continue;
 
        // get global id for the element
        UInt glb = mesh.elementToLinearized(facet_elems(i)[j]);
        
        // check if the vertex corresponding to the element has been inserted
        check_iterator uit = inserted.find(glb);
        if (uit == inserted.end()) {
          u.push_back(boost::add_vertex(g_));
          index[u[j]] = glb;
          inserted[glb] = u[j];
        } else
          u.push_back(inserted[glb]);
      }
      // if facet connected to two elements, add edge
      if (u.size() == 2)
        boost::add_edge(u[0],u[1],g_);
    }
  }
  
  // call connected components algorithm from boost
  connectedComponents(mesh);
}

//MeshGraph::MeshGraph(Mesh& mesh, Mesh& facets) : components_("components","mesh_graph") {
//  
//  // get mesh
//  UInt d = mesh.getSpatialDimension();
//  
//  mesh.initElementTypeMapArray(components_,1,d,false,_ek_regular,true);
//
//  // create graph
//  g_ = graph_type();
//
//  // obtain a property map for the vertex_index property
//  boost::property_map<graph_type, boost::vertex_index1_t>::type index = boost::get(boost::vertex_index1_t(), g_);
//
//
//  std::map<UInt, vertex_type> vertex_map;
//  
//  UInt total_els = 0;
//  for (Mesh::type_iterator it = mesh.firstType(d); it != mesh.lastType(d); ++it) {
////    total_els += mesh.getNbElement(*it);
//
//    UInt els = mesh.getNbElement(*it);
//    Element el(*it);
//    
//    for (UInt iel = 0; iel<els; ++iel) {
//      
//      el.element = iel;
//      vertex_type u = boost::add_vertex(g_);
//      UInt glb = mesh.elementToLinearized(el);
//      index[u] = glb;
//      vertex_map[glb] = u;
////      boost::put(index, u, glb);
//    }
//  }
//  
//  // loop over facets
//  for (auto it_f = facets.firstType(d - 1); it_f != facets.lastType(d - 1); ++it_f) {
//    
//    UInt nb_facet = facets.getNbElement(*it_f);
//    
//    /// get the element connected to a subelement
//    const Array< std::vector<Element> > &facet_elems = facets.getElementToSubelement(*it_f);
//    
//    AKANTU_DEBUG_ASSERT(facet_elems.getSize() == nb_facet, "Number of facets does not match array dimentions");
//    
//    // loop over facets to construct mesh dual graph
//    for (UInt i=0; i<nb_facet; ++i) {
//      
//      // skip if facet is on the boundary
//      if (facet_elems(i)[0].ghost_type == _ghost || facet_elems(i)[1].ghost_type == _ghost || facet_elems(i)[1] == ElementNull)
//        continue;
//      
//      cout<<"facet_elems(i)[0].element -> "<<(facet_elems(i)[0].element)<<endl;
//      
//      cout<<"mesh.elementToLinearized(facet_elems(i)[0]) -> "<<mesh.elementToLinearized(facet_elems(i)[0])<<endl;
//      
//      // add edge
////      boost::add_edge(facet_elems(i)[0].element,
////                      facet_elems(i)[1].element,
////                      g_);
////      boost::add_edge(mesh.elementToLinearized(facet_elems(i)[0]),
////                      mesh.elementToLinearized(facet_elems(i)[1]),
////                      g_);
//      boost::add_edge(vertex_map[mesh.elementToLinearized(facet_elems(i)[0])],
//                      vertex_map[mesh.elementToLinearized(facet_elems(i)[1])],
//                      g_);
//    }
//  }
//
//  connectedComponents(mesh);
//}



//template <class filter_type>
//MeshGraph::MeshGraph(Mesh& mesh, Mesh& facets, filter_type& filter) : components_("components","mesh_graph") {
//  
//  // get mesh
//  UInt d = mesh.getSpatialDimension();
//  
//  mesh.initElementTypeMapArray(components_,1,d,false,_ek_regular,true);
//  
//  UInt total_els = 0;
//  for (Mesh::type_iterator it = mesh.firstType(d,_not_ghost,_ek_not_defined); it != mesh.lastType(d,_not_ghost,_ek_not_defined); ++it) {
//    auto el_filter = filter(*it);
//    for (size_t i=0; i<el_filter.getSize(); ++i)
//      total_els += el_filter(i);
//  }
//  
//  // create graph with the total number of elements
//  g_ = graph_type(total_els);
//  
//  // loop over facets
//  for (auto it_f = facets.firstType(d - 1); it_f != facets.lastType(d - 1); ++it_f) {
//    
//    UInt nb_facet = facets.getNbElement(*it_f);
//    
//    /// get the element connected to a subelement
//    const Array< std::vector<Element> > &facet_elems = facets.getElementToSubelement(*it_f);
//    
//    AKANTU_DEBUG_ASSERT(facet_elems.getSize() == nb_facet, "Number of facets does not match array dimentions");
//    
//    // loop over facets to construct mesh dual graph
//    for (UInt i=0; i<nb_facet; ++i) {
//      
//      bool skip = false;
//      // skip with element filter
//      for (size_t j = 0; j<2; ++j) {
//        const Element &el = facet_elems(i)[j];
//        if (el.kind == _ek_cohesive) {
//          if (filter(el.element))
//            skip = true;
//        }
//        // skip if facet is on the boundary
//        else if (el.ghost_type == _ghost || el == ElementNull)
//          skip = true;
//      }
//      
//      if (skip) continue;
//      
////      // add connectivity
////      boost::add_edge(mesh.elementToLinearized(facet_elems(i)[0]),
////                      mesh.elementToLinearized(facet_elems(i)[1]),
////                      g_);
//    }
//  }
//  
//  connectedComponents(mesh);
//  
//}




__END_AKANTU__
