/**
 * @file   mesh.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Jun 18 11:47:19 2010
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

#include "aka_config.hh"

#ifdef AKANTU_USE_CPPARRAY
#include <array/expr.hpp>
#endif

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "mesh_io.hh"
#include "element_class.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

const Element ElementNull(_not_defined, 0);

/* -------------------------------------------------------------------------- */
void Element::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  stream << space << "Element [" << type << ", " << element << ", " << ghost_type << "]";
}


/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
           const ID id,
           const MemoryID & memory_id) :
  Memory(memory_id), id(id), nodes_global_ids(NULL), nodes_type(NULL),
  created_nodes(true),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  ghost_types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  mesh_data("mesh_data", id, memory_id),
  boundaries(*this, "boundaries", id, memory_id),
  facet_to_double("facet_to_double", id),
  subfacets_to_subsubfacet_double("subfacets_to_subsubfacet_double", id),
  facets_to_subfacet_double("facets_to_subfacet_double", id),
  elements_to_subfacet_double("elements_to_subfacet_double", id) {
  AKANTU_DEBUG_IN();

  this->nodes = &(alloc<Real>(this->id + ":coordinates", 0, this->spatial_dimension));

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
  types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  ghost_types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  mesh_data("mesh_data", id, memory_id),
  boundaries(*this, "boundaries", id, memory_id),
  facet_to_double("facet_to_double", id),
  subfacets_to_subsubfacet_double("subfacets_to_subsubfacet_double", id),
  facets_to_subfacet_double("facets_to_subfacet_double", id),
  elements_to_subfacet_double("elements_to_subfacet_double", id) {
  AKANTU_DEBUG_IN();

  this->nodes = &(getArray<Real>(nodes_id));
  nb_global_nodes = nodes->getSize();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
           Array<Real> & nodes,
           const ID & id,
           const MemoryID & memory_id) :
  Memory(memory_id), id(id), nodes_global_ids(NULL), nodes_type(NULL),
  created_nodes(false),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Array<UInt>(_max_element_type + 1, 1)),
  ghost_types_offsets(Array<UInt>(_max_element_type + 1, 1)),
  mesh_data("mesh_data", id, memory_id),
  boundaries(*this, "boundaries", id, memory_id),
  facet_to_double("facet_to_double", id),
  subfacets_to_subsubfacet_double("subfacets_to_subsubfacet_double", id),
  facets_to_subfacet_double("facets_to_subfacet_double", id),
  elements_to_subfacet_double("elements_to_subfacet_double", id) {
  AKANTU_DEBUG_IN();

  this->nodes = &(nodes);
  nb_global_nodes = nodes.getSize();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh & Mesh::initMeshFacets(const ID & id) {
  AKANTU_DEBUG_IN();

  if (!mesh_facets) {
    mesh_facets = new Mesh(spatial_dimension,
                           *(this->nodes),
                           this->id+":"+id,
                           getMemoryID());

    mesh_facets->mesh_parent = this;
    mesh_facets->is_mesh_facets = true;
  }

  AKANTU_DEBUG_OUT();
  return *mesh_facets;
}

/* -------------------------------------------------------------------------- */
void Mesh::init() {
  nodes_type = NULL;
  mesh_facets = NULL;
  is_mesh_facets = false;
  mesh_parent = NULL;
  //  computeBoundingBox();
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {
  AKANTU_DEBUG_IN();

  delete mesh_facets;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::read (const std::string & filename, const MeshIOType & mesh_io_type) {
  MeshIO mesh_io;
  mesh_io.read(filename, *this, mesh_io_type);
}

/* -------------------------------------------------------------------------- */
void Mesh::write(const std::string & filename, const MeshIOType & mesh_io_type) {
  MeshIO mesh_io;
  mesh_io.write(filename, *this, mesh_io_type);
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
    //    if(!isPureGhostNode(i))
    for (UInt k = 0; k < spatial_dimension; ++k) {
      local_lower_bounds[k] = std::min(local_lower_bounds[k], (*nodes)(i, k));
      local_upper_bounds[k] = std::max(local_upper_bounds[k], (*nodes)(i, k));
    }
  }

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  Real reduce_bounds[2 * spatial_dimension];
  for (UInt k = 0; k < spatial_dimension; ++k) {
    reduce_bounds[2*k    ] =   local_lower_bounds[k];
    reduce_bounds[2*k + 1] = - local_upper_bounds[k];
  }

  comm.allReduce(reduce_bounds, 2 * spatial_dimension, _so_min);

  for (UInt k = 0; k < spatial_dimension; ++k) {
    lower_bounds[k] =   reduce_bounds[2*k];
    upper_bounds[k] = - reduce_bounds[2*k + 1];
  }

  for (UInt k = 0; k < spatial_dimension; ++k)
    size[k] = upper_bounds[k] - lower_bounds[k];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Mesh::initByElementTypeArray(ByElementTypeArray<T> & vect,
                                  UInt nb_component,
                                  UInt dim,
                                  const bool & flag_nb_node_per_elem_multiply,
                                  ElementKind element_kind,
                                  bool size_to_nb_element) const {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    this->initByElementTypeArray(vect, nb_component, dim, gt,
                                 flag_nb_node_per_elem_multiply,
                                 element_kind, size_to_nb_element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Mesh::initByElementTypeArray(ByElementTypeArray<T> & vect,
                                  UInt nb_component,
                                  UInt dim,
                                  GhostType gt,
                                  const bool & flag_nb_node_per_elem_multiply,
                                  ElementKind element_kind,
                                  bool size_to_nb_element) const {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it  = firstType(dim, gt, element_kind);
  Mesh::type_iterator end = lastType(dim, gt, element_kind);
  for(; it != end; ++it) {
    ElementType type = *it;
    if (flag_nb_node_per_elem_multiply) nb_component *= Mesh::getNbNodesPerElement(*it);
    UInt size = 0;
    if (size_to_nb_element) size = this->getNbElement(type, gt);
    vect.alloc(size, nb_component, type, gt);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::initNormals() {
  initByElementTypeArray(normals, spatial_dimension, spatial_dimension, false, _ek_not_defined);
}

/* -------------------------------------------------------------------------- */
void Mesh::initFacetToDouble() {
  for (UInt sp = 0; sp < spatial_dimension; ++sp) {
    initByElementTypeArray(facet_to_double, 1, sp);
    initByElementTypeArray(facets_to_subfacet_double, 1, sp);
    initByElementTypeArray(elements_to_subfacet_double, 1, sp);
    initByElementTypeArray(subfacets_to_subsubfacet_double, 1, sp);
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::getGlobalConnectivity(Array<UInt> & global_connectivity,
				 ElementType type,
				 GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<UInt> & local_conn = connectivities(type, ghost_type);
  if (!nodes_global_ids)
    nodes_global_ids = mesh_parent->nodes_global_ids;

  UInt * local_c = local_conn.storage();
  UInt * global_c = global_connectivity.storage();

  UInt nb_terms = local_conn.getSize() * local_conn.getNbComponent();

  for (UInt i = 0; i < nb_terms; ++i, ++local_c, ++global_c)
    *global_c = (*nodes_global_ids)(*local_c);

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
template void Mesh::initByElementTypeArray<Real>(ByElementTypeArray<Real> & vect,
                                                 UInt nb_component,
                                                 UInt dim,
                                                 const bool & flag_nb_elem_multiply,
                                                 ElementKind element_kind,
                                                 bool size_to_nb_element) const;

template void Mesh::initByElementTypeArray<Int>(ByElementTypeArray<Int> & vect,
                                                UInt nb_component,
                                                UInt dim,
                                                const bool & flag_nb_elem_multiply,
                                                ElementKind element_kind,
                                                bool size_to_nb_element) const;

template void Mesh::initByElementTypeArray<UInt>(ByElementTypeArray<UInt> & vect,
                                                 UInt nb_component,
                                                 UInt dim,
                                                 const bool & flag_nb_elem_multiply,
                                                 ElementKind element_kind,
                                                 bool size_to_nb_element) const;

template void Mesh::initByElementTypeArray<bool>(ByElementTypeArray<bool> & vect,
                                                 UInt nb_component,
                                                 UInt dim,
                                                 const bool & flag_nb_elem_multiply,
                                                 ElementKind element_kind,
                                                 bool size_to_nb_element) const;


__END_AKANTU__
