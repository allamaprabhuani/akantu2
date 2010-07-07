/**
 * @file   mesh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 12:02:26 2010
 *
 * @brief  class handling meshes
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <sstream>

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

Mesh::Mesh(unsigned int spatial_dimension, MemoryID memory_id) :
  Memory(memory_id) {

  this->spatial_dimension = spatial_dimension;
  nodes = NULL;
}

/* -------------------------------------------------------------------------- */

Vector<int> & Mesh::createConnectivity(ElementType type, unsigned int nb_element) {

  unsigned int nb_nodes_per_element;

  switch(type) {
  case _triangle_1   : {
    ElementClass<_triangle_1>   elem_class;
    nb_nodes_per_element = elem_class.getNbNodesPerElement();
    break; }
  case _triangle_2   : {
    ElementClass<_triangle_2>   elem_class;
    nb_nodes_per_element = elem_class.getNbNodesPerElement();
    break; }
  case _tetrahedra_1 : {
    ElementClass<_tetrahedra_1> elem_class;
    nb_nodes_per_element = elem_class.getNbNodesPerElement();
    break; }
  case _tetrahedra_2 : {
    ElementClass<_tetrahedra_2> elem_class;
    nb_nodes_per_element = elem_class.getNbNodesPerElement();
    break; }
  case _not_defined :
  case _max_element_type : {
    nb_nodes_per_element = 0;
    AKANTU_DEBUG_ERROR("Cannot create a conectivity vector of type : " << type);
    break;
  }
  }

  std::stringstream sstr;
  sstr << id << ":connectivity:" << type;

  Vector<int> & connectivity = malloc<int>(sstr.str(),
					   nb_element,
					   nb_nodes_per_element);

  connectivities[type] = &connectivity;
  return connectivity;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
