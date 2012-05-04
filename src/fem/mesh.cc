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
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

const Element ElementNull(_not_defined, 0);

/* -------------------------------------------------------------------------- */
void Element::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  stream << space << "Element [" << type << ", " << element << "]";
}


/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   const ID & id,
	   const MemoryID & memory_id) :
  Memory(memory_id), id(id), nodes_global_ids(NULL), nodes_type(NULL),
  created_nodes(true),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Vector<UInt>((UInt) _max_element_type + 1, 1)),
  ghost_types_offsets(Vector<UInt>((UInt) _max_element_type + 1, 1)),
  nb_surfaces(0),
  surface_id("surface_id", id),
  element_to_subelement("element_to_subelement", id),
  subelement_to_element("subelement_to_element", id),
  uint_data("by_element_uint_data", id) {
  AKANTU_DEBUG_IN();

  std::stringstream sstr;
  sstr << id << ":coordinates";
  this->nodes = &(alloc<Real>(sstr.str(), 0, this->spatial_dimension));

  nb_global_nodes = 0;

  init();

  std::fill_n(lower_bounds, 3, 0.);
  std::fill_n(upper_bounds, 3, 0.);

  std::fill_n(size, 3, 0.);

  std::fill_n(local_lower_bounds, 3, 0.);
  std::fill_n(local_upper_bounds, 3, 0.);

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   const ID & nodes_id,
	   const ID & id,
	   const MemoryID & memory_id) :
  Memory(memory_id), id(id), nodes_global_ids(NULL), nodes_type(NULL),
  created_nodes(false),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Vector<UInt>((UInt) _max_element_type + 1, 1)),
  ghost_types_offsets(Vector<UInt>((UInt) _max_element_type + 1, 1)),
  nb_surfaces(0),
  surface_id("surface_id", id),
  element_to_subelement("element_to_subelement", id),
  subelement_to_element("subelement_to_element", id),
  uint_data("by_element_uint_data", id) {
  AKANTU_DEBUG_IN();

  this->nodes = &(getVector<Real>(nodes_id));
  nb_global_nodes = nodes->getSize();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   Vector<Real> & nodes,
	   const ID & id,
	   const MemoryID & memory_id) :
  Memory(memory_id), id(id), nodes_global_ids(NULL), nodes_type(NULL),
  created_nodes(false),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Vector<UInt>(_max_element_type + 1, 1)),
  ghost_types_offsets(Vector<UInt>(_max_element_type + 1, 1)),
  nb_surfaces(0),
  surface_id("surface_id", id),
  element_to_subelement("element_to_subelement", id),
  subelement_to_element("subelement_to_element", id),
  uint_data("by_element_uint_data", id) {
  AKANTU_DEBUG_IN();

  this->nodes = &(nodes);
  nb_global_nodes = nodes.getSize();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::init() {
  //  this->types_offsets.resize(_max_element_type);

  nodes_type = NULL;
  computeBoundingBox();
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it  = firstType(0, gt);
    Mesh::type_iterator end = lastType(0, gt);
    for(; it != end; ++it) {
      UIntDataMap & map = uint_data(*it, gt);
      UIntDataMap::iterator dit;
      for (dit = map.begin(); dit != map.end(); ++dit) {
	if(dit->second) delete dit->second;
      }
      map.clear();
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

  stream << space << " + connectivities [" << std::endl;
  connectivities.printself(stream, indent+2);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void Mesh::computeBoundingBox(){
  AKANTU_DEBUG_IN();
  for (UInt k = 0; k < spatial_dimension; ++k) {
    local_lower_bounds[k] =   std::numeric_limits<double>::max();
    local_upper_bounds[k] = - std::numeric_limits<double>::max();
  }

  for (UInt i = 0; i < nodes->getSize(); ++i) {
    if(!isPureGhostNode(i))
      for (UInt k = 0; k < spatial_dimension; ++k) {
        local_lower_bounds[k] = std::min(local_lower_bounds[k], (*nodes)(i, k));
        local_upper_bounds[k] = std::max(local_upper_bounds[k], (*nodes)(i, k));
      }
  }

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();

  Real reduce_bounds[2 * spatial_dimension];
  for (UInt k = 0; k < spatial_dimension; ++k) {
    reduce_bounds[2*k    ] =   local_lower_bounds[k];
    reduce_bounds[2*k + 1] = - local_upper_bounds[k];
  }

  comm->allReduce(reduce_bounds, 2 * spatial_dimension, _so_min);

  for (UInt k = 0; k < spatial_dimension; ++k) {
    lower_bounds[k] =   reduce_bounds[2*k];
    upper_bounds[k] = - reduce_bounds[2*k + 1];
  }

  for (UInt k = 0; k < spatial_dimension; ++k)
    size[k] = upper_bounds[k] - lower_bounds[k];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::setSurfaceIDsFromIntData(const std::string & data_name) {

  std::set<Surface> surface_ids;

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it  = firstType(spatial_dimension - 1, gt);
    Mesh::type_iterator end = lastType(spatial_dimension - 1, gt);
    for(; it != end; ++it) {
      UIntDataMap & map = uint_data(*it, gt);
      UIntDataMap::iterator it_data = map.find(data_name);
      AKANTU_DEBUG_ASSERT(it_data != map.end(),
			  "No data named " << data_name
			  << " present in the mesh " << id
			  << " for the element type " << *it);
      AKANTU_DEBUG_ASSERT(!surface_id.exists(*it, gt),
			  "Surface id for type (" << gt << ":" << *it
			  << ") already set to the vector " << surface_id(*it, gt).getID());

      surface_id.setVector(*it, gt, *it_data->second);

      for (UInt s = 0; s < it_data->second->getSize(); ++s) {
	surface_ids.insert((*it_data->second)(s));
      }
    }
  }

  nb_surfaces = surface_ids.size();
}


/* -------------------------------------------------------------------------- */
template<typename T>
void Mesh::initByElementTypeVector(ByElementTypeVector<T> & vect,
				   UInt nb_component,
				   UInt dim,
				   const bool & flag_nb_node_per_elem_multiply,
				   ElementKind element_kind) const {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it  = firstType(dim, gt, element_kind);
    Mesh::type_iterator end = lastType(dim, gt, element_kind);
    for(; it != end; ++it) {
      ElementType type = *it;
      if (flag_nb_node_per_elem_multiply) nb_component *= Mesh::getNbNodesPerElement(*it);
      vect.alloc(0, nb_component, type);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template void Mesh::initByElementTypeVector<Real>(ByElementTypeVector<Real> & vect,
						  UInt nb_component,
						  UInt dim,
						  const bool & flag_nb_elem_multiply,
						  ElementKind element_kind) const;

template void Mesh::initByElementTypeVector<Int>(ByElementTypeVector<Int> & vect,
						 UInt nb_component,
						 UInt dim,
						 const bool & flag_nb_elem_multiply,
						 ElementKind element_kind) const;

template void Mesh::initByElementTypeVector<UInt>(ByElementTypeVector<UInt> & vect,
						  UInt nb_component,
						  UInt dim,
						  const bool & flag_nb_elem_multiply,
						  ElementKind element_kind) const;

__END_AKANTU__
