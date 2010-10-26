/**
 * @file   mesh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 12:02:26 2010
 *
 * @brief  class handling meshes
 *
 * @section LICENSE
 *
 * \<insert license here\>
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
void Element::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  stream << space << "Element [" << type << ", " << element << "]";
}


/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   const MeshID & id,
	   const MemoryID & memory_id) :
  Memory(memory_id), id(id), nodes_global_ids(NULL),
  created_nodes(true), spatial_dimension(spatial_dimension),
  internal_facets_mesh(NULL),
  types_offsets(Vector<UInt>(_max_element_type + 1, 1)),
  ghost_types_offsets(Vector<UInt>(_max_element_type + 1, 1)) {
  AKANTU_DEBUG_IN();

  initConnectivities();

  std::stringstream sstr;
  sstr << id << ":coordinates";
  this->nodes = &(alloc<Real>(sstr.str(), 0, this->spatial_dimension));

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   const VectorID & nodes_id,
	   const MeshID & id,
	   const MemoryID & memory_id) :
  Memory(memory_id), id(id), nodes_global_ids(NULL),
  created_nodes(false), spatial_dimension(spatial_dimension),
  internal_facets_mesh(NULL),
  types_offsets(Vector<UInt>(_max_element_type + 1, 1)),
  ghost_types_offsets(Vector<UInt>(_max_element_type + 1, 1)) {
  AKANTU_DEBUG_IN();

  initConnectivities();

  this->nodes = &(getVector<Real>(nodes_id));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   Vector<Real> & nodes,
	   const MeshID & id,
	   const MemoryID & memory_id) :
  Memory(memory_id), id(id), created_nodes(false), spatial_dimension(spatial_dimension),
  internal_facets_mesh(NULL),
  types_offsets(Vector<UInt>(_max_element_type + 1, 1)),
  ghost_types_offsets(Vector<UInt>(_max_element_type + 1, 1)) {
  AKANTU_DEBUG_IN();

  initConnectivities();

  this->nodes = &(nodes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::initConnectivities() {
  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    connectivities[t] = NULL;
    ghost_connectivities[t] = NULL;
    surface_id[t] = NULL;
  }
  this->types_offsets.resize(_max_element_type);
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Vector<Real> & Mesh::createNormals(ElementType type) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ERROR("TOBEIMPLEMENTED");
  AKANTU_DEBUG_OUT();
  return *normals[type];
}

/* -------------------------------------------------------------------------- */
void Mesh::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Mesh [" << std::endl;
  stream << space << " + id                : " << this->id << std::endl;
  stream << space << " + spatial dimension : " << this->spatial_dimension << std::endl;
  stream << space << " + nodes [" << std::endl;
  nodes->printself(stream, indent+2);
  stream << space << " ]" << std::endl;

  ConnectivityTypeList::const_iterator it;
  for(it = type_set.begin(); it != type_set.end(); ++it) {
    stream << space << " + connectivities ("<< *it <<") [" << std::endl;
    (connectivities[*it])->printself(stream, indent+2);
    stream << space << " ]" << std::endl;
  }
  for(it = ghost_type_set.begin(); it != ghost_type_set.end(); ++it) {
    stream << space << " + ghost_connectivities ("<< *it <<") [" << std::endl;
    (ghost_connectivities[*it])->printself(stream, indent+2);
    stream << space << " ]" << std::endl;
  }
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
