/**
 * @file   mesh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 12:02:26 2010
 *
 * @brief  class handling meshes
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <sstream>

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   const MeshID & id,
	   const MemoryID & memory_id) :
  Memory(memory_id), id(id), created_nodes(true), spatial_dimension(spatial_dimension) {
  AKANTU_DEBUG_IN();

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    connectivities[t] = NULL;
  }

  std::stringstream sstr;
  sstr << id << ":coordinates";
  nodes = &(alloc<double>(sstr.str(), 0, this->spatial_dimension));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {
  AKANTU_DEBUG_IN();
  if(created_nodes) {
    AKANTU_DEBUG(dblAccessory, "Deleting nodes vector");
    dealloc(nodes->getID());
  }
  nodes = NULL;

  ConnectivityTypeList::const_iterator it;
  for(it = type_set.begin();
      it != type_set.end();
      ++it) {
    AKANTU_DEBUG(dblAccessory, "Deleting connectivity vector of type " << *it);
    dealloc(connectivities[*it]->getID());
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Vector<UInt> & Mesh::createConnectivity(ElementType type, UInt nb_element) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes_per_element;

  AKANTU_DEBUG_ASSERT(connectivities[type] == NULL,
		      "The connectivity vector for the type "
		      << type << "already exist");

#define GET_NB_NODES_PER_ELEM(type)					\
  nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement()

  switch(type) {
  case _line_1       : { GET_NB_NODES_PER_ELEM(_line_1      ); break; }
  case _line_2       : { GET_NB_NODES_PER_ELEM(_line_2      ); break; }
  case _triangle_1   : { GET_NB_NODES_PER_ELEM(_triangle_1  ); break; }
  case _triangle_2   : { GET_NB_NODES_PER_ELEM(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_NB_NODES_PER_ELEM(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_NB_NODES_PER_ELEM(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }
#undef GET_NB_NODES_PER_ELEM

  std::stringstream sstr;
  sstr << id << ":connectivity:" << type;

  connectivities[type] = &(alloc<UInt>(sstr.str(),
				     nb_element,
				     nb_nodes_per_element));

  type_set.insert(type);

  AKANTU_DEBUG_OUT();
  return *connectivities[type];
}

/* -------------------------------------------------------------------------- */
void Mesh::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Mesh [" << std::endl;
  stream << space << " + id             : " << this->id << std::endl;
  stream << space << " + spatial dim    : " << this->spatial_dimension << std::endl;
  stream << space << " + nodes [" << std::endl;
  nodes->printself(stream, indent+2);
  stream << space << " ]" << std::endl;

  ConnectivityTypeList::const_iterator it;
  for(it = type_set.begin(); it != type_set.end(); ++it) {
    stream << space << " + connectivities ("<< *it <<") [" << std::endl;
    (connectivities[*it])->printself(stream, indent+2);
    stream << space << " ]" << std::endl;
  }
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
