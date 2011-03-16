/**
 * @file   mesh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 12:02:26 2010
 *
 * @brief  class handling meshes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
  ghost_types_offsets(Vector<UInt>(_max_element_type + 1, 1)),
  nb_surfaces(0) {
  AKANTU_DEBUG_IN();

  std::stringstream sstr;
  sstr << id << ":coordinates";
  this->nodes = &(alloc<Real>(sstr.str(), 0, this->spatial_dimension));

  nb_global_nodes = 0;

  init();

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

  this->nodes = &(getVector<Real>(nodes_id));
  nb_global_nodes = nodes->getSize();

  init();

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

  this->nodes = &(nodes);
  nb_global_nodes = nodes.getSize();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::init() {
  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    connectivities[t] = NULL;
    ghost_connectivities[t] = NULL;
    surface_id[t] = NULL;

    uint_data[t].clear();
    ghost_uint_data[t].clear();
  }

  this->types_offsets.resize(_max_element_type);

  nodes_type = NULL;
  computeBoundingBox();
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {
  AKANTU_DEBUG_IN();

  for (UInt t = _not_defined;  t < _max_element_type; ++t) {
    if(uint_data[t].size() > 0) {
      UIntDataMap::iterator it;
      for (it = uint_data[t].begin(); it != uint_data[t].end(); ++it) {
	if(it->second) delete it->second;
      }
      uint_data[t].clear();
    }
    if(ghost_uint_data[t].size() > 0) {
      UIntDataMap::iterator it;
      for (it = ghost_uint_data[t].begin(); it != ghost_uint_data[t].end(); ++it) {
	if(it->second) delete it->second;
      }
      ghost_uint_data[t].clear();
    }
  }

  AKANTU_DEBUG_OUT();
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

void Mesh::computeBoundingBox(){
  AKANTU_DEBUG_IN();
  UInt dim = spatial_dimension;

  for (UInt k = 0; k < dim; ++k) {
    xmin[k] = std::numeric_limits<double>::max();
    xmax[k] = std::numeric_limits<double>::min();
  }

  Real * coords = nodes->values;
  for (UInt i = 0; i < nodes->getSize(); ++i) {
    for (UInt k = 0; k < dim; ++k) {
      xmin[k] = std::min(xmin[k],coords[dim*i+k]);
      xmax[k] = std::max(xmax[k],coords[dim*i+k]);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::initByElementTypeRealVector(ByElementTypeReal & vect,
				       UInt nb_component,
				       UInt dim,
				       const std::string & obj_id,
				       const std::string & vect_id,
				       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for(UInt t = _not_defined; t < _max_element_type; ++t)
    vect[t] = NULL;

  std::string ghost_id = "";

  if (ghost_type == _ghost) {
    ghost_id = "ghost_";
  }

  const Mesh::ConnectivityTypeList & type_list = getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(dim > 0 && Mesh::getSpatialDimension(*it) != dim) continue;
    std::stringstream sstr; sstr << obj_id << ":" << ghost_id << vect_id << ":" << *it;
    if (vect[*it] == NULL){
      vect[*it] = &(alloc<Real>(sstr.str(), 0,
				nb_component, REAL_INIT_VALUE));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::initByElementTypeUIntVector(ByElementTypeUInt & vect,
				       UInt nb_component,
				       UInt dim,
				       const std::string & obj_id,
				       const std::string & vect_id,
				       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for(UInt t = _not_defined; t < _max_element_type; ++t)
    vect[t] = NULL;

  std::string ghost_id = "";

  if (ghost_type == _ghost) {
    ghost_id = "ghost_";
  }

  const Mesh::ConnectivityTypeList & type_list = getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(dim > 0 && Mesh::getSpatialDimension(*it) != dim) continue;
    std::stringstream sstr; sstr << obj_id << ":" << ghost_id << vect_id << ":" << *it;
    if (vect[*it] == NULL){
      vect[*it] = &(alloc<UInt>(sstr.str(), 0,
				nb_component, UINT_INIT_VALUE));
    }
  }

  AKANTU_DEBUG_OUT();
}




__END_AKANTU__
