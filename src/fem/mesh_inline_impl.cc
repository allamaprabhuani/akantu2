/**
 * @file   mesh_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Jul 15 00:41:12 2010
 *
 * @brief  Implementation of the inline functions of the mesh class
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

__END_AKANTU__
#if defined(AKANTU_COHESIVE_ELEMENT)
#  include "cohesive_element.hh"
#endif
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline RemovedNodesEvent::RemovedNodesEvent(const Mesh & mesh) :
  new_numbering(mesh.getNbNodes(), 1, "new_numbering") {
}

/* -------------------------------------------------------------------------- */
inline RemovedElementsEvent::RemovedElementsEvent(const Mesh & mesh) :
  new_numbering("new_numbering", mesh.getID()) {
}

/* -------------------------------------------------------------------------- */
template <>
inline void Mesh::sendEvent<RemovedElementsEvent>(RemovedElementsEvent & event) {
  connectivities.onElementsRemoved(event.getNewNumbering());
  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <>
inline void Mesh::sendEvent<RemovedNodesEvent>(RemovedNodesEvent & event) {
  if(created_nodes)    removeNodesFromVector(*nodes           , event.getNewNumbering());
  if(nodes_global_ids) removeNodesFromVector(*nodes_global_ids, event.getNewNumbering());
  if(nodes_type)       removeNodesFromVector(*nodes_type      , event.getNewNumbering());

  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Mesh::removeNodesFromVector(Vector<T> & vect, const Vector<UInt> & new_numbering) {
  Vector<T> tmp(vect.getSize(), vect.getNbComponent());
  UInt nb_component = vect.getNbComponent();
  UInt new_nb_nodes = 0;
  for (UInt i = 0; i < new_numbering.getSize(); ++i) {
    UInt new_i = new_numbering(i);
    if(new_i != UInt(-1)) {
      memcpy(tmp.storage() + new_i * nb_component,
	     vect.storage() + i * nb_component,
	     nb_component * sizeof(T));
      ++new_nb_nodes;
    }
  }

  tmp.resize(new_nb_nodes);
  vect.copy(tmp);
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::elementToLinearized(const Element & elem) const {
  AKANTU_DEBUG_ASSERT(elem.type < _max_element_type &&
		      elem.element < types_offsets.values[elem.type+1],
		      "The element " << elem
		      << "does not exists in the mesh " << id);

  return types_offsets.values[elem.type] + elem.element;
}

/* -------------------------------------------------------------------------- */
inline Element Mesh::linearizedToElement (UInt linearized_element) const {

  UInt t;

  for (t = _not_defined;
       t != _max_element_type && linearized_element >= types_offsets(t);
       ++t);

  AKANTU_DEBUG_ASSERT(linearized_element < types_offsets(t),
   		      "The linearized element " << linearized_element
   		      << "does not exists in the mesh " << id);

  --t;
  ElementType type = ElementType(t);
  return Element(type,
		 linearized_element - types_offsets.values[t],
		 _not_ghost,
		 getKind(type));
}

/* -------------------------------------------------------------------------- */
inline void Mesh::updateTypesOffsets(const GhostType & ghost_type) {
  types_offsets.clear();
  type_iterator it   = firstType(0, ghost_type, _ek_not_defined);
  type_iterator last = lastType(0, ghost_type, _ek_not_defined);

  for (; it != last; ++it)
    types_offsets(*it) = connectivities(*it, ghost_type).getSize();

  for (UInt t = _not_defined + 1; t < _max_element_type; ++t)
    types_offsets(t) += types_offsets(t - 1);
  for (UInt t = _max_element_type; t > _not_defined; --t)
    types_offsets(t) = types_offsets(t - 1);
  types_offsets(0) = 0;
}

/* -------------------------------------------------------------------------- */
inline const Mesh::ConnectivityTypeList & Mesh::getConnectivityTypeList(const GhostType & ghost_type) const {
  if (ghost_type == _not_ghost)
    return type_set;
  else
    return ghost_type_set;
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getNodesGlobalIdsPointer() {
  AKANTU_DEBUG_IN();
  if(nodes_global_ids == NULL) {
    std::stringstream sstr; sstr << id << ":nodes_global_ids";
    nodes_global_ids = &(alloc<UInt>(sstr.str(), nodes->getSize(), 1));
  }
  AKANTU_DEBUG_OUT();
  return nodes_global_ids;
}

/* -------------------------------------------------------------------------- */
inline Vector<Int> * Mesh::getNodesTypePointer() {
  AKANTU_DEBUG_IN();
  if(nodes_type == NULL) {
    std::stringstream sstr; sstr << id << ":nodes_type";
    nodes_type = &(alloc<Int>(sstr.str(), nodes->getSize(), 1, -1));
  }
  AKANTU_DEBUG_OUT();
  return nodes_type;
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getConnectivityPointer(const ElementType & type,
						   const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Vector<UInt> * tmp;
  if(!connectivities.exists(type, ghost_type)) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    tmp = &(connectivities.alloc(0, nb_nodes_per_element,
				 type, ghost_type));

    AKANTU_DEBUG_INFO("The connectivity vector for the type "
		      << type << " created");

    if (ghost_type == _not_ghost) type_set.insert(type);
    else ghost_type_set.insert(type);

    updateTypesOffsets(ghost_type);
  } else {
    tmp = &connectivities(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
  return tmp;
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getSurfaceIDPointer(const ElementType & type,
						const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Vector<UInt> * tmp;
  if(!surface_id.exists(type, ghost_type)) {
    tmp = &(surface_id.alloc(0, 1, type, ghost_type));
    AKANTU_DEBUG_INFO("The surface id vector for the type "
		      << type << " created");
  } else {
    tmp = &(surface_id(type, ghost_type));
  }

  AKANTU_DEBUG_OUT();
  return tmp;
}

/* -------------------------------------------------------------------------- */
inline Vector< std::vector<Element> > * Mesh::getElementToSubelementPointer(const ElementType & type,
									    const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Vector< std::vector<Element> > * tmp;
  if(!element_to_subelement.exists(type, ghost_type)) {
    tmp = &(element_to_subelement.alloc(0, 1, type, ghost_type));

    AKANTU_DEBUG_INFO("The element_to_subelement vector for the type "
		      << type << " created");
  } else {
    tmp = &(element_to_subelement(type, ghost_type));
  }

  AKANTU_DEBUG_OUT();
  return tmp;
}

/* -------------------------------------------------------------------------- */
inline Vector<Element > * Mesh::getSubelementToElementPointer(const ElementType & type,
							      const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Vector<Element > * tmp;
  if(!subelement_to_element.exists(type, ghost_type)) {

    UInt nb_facets = getNbFacetsPerElement(type);

    tmp = &(subelement_to_element.alloc(0, nb_facets, type, ghost_type));

    AKANTU_DEBUG_INFO("The subelement_to_element vector for the type "
		      << type << " created");
  } else {
    tmp = &(subelement_to_element(type, ghost_type));
  }

  AKANTU_DEBUG_OUT();
  return tmp;
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getUIntDataPointer(const ElementType & el_type,
					       const std::string & data_name,
					       const GhostType & ghost_type) {
  //  AKANTU_DEBUG_IN();

  Vector<UInt> * data;
  // if(!uint_data.exists(el_type, ghost_type)){
  //   uint_data(UIntDataMap(), el_type, ghost_type);
  // }
  UIntDataMap & map = uint_data(el_type, ghost_type);
  UIntDataMap::iterator it = map.find(data_name);
  if(it == map.end()) {
    data = new Vector<UInt>(0, 1, data_name);
    map[data_name] = data;
  } else {
    data = it->second;
  }

  //  AKANTU_DEBUG_OUT();
  return data;
}

/* -------------------------------------------------------------------------- */
inline const Vector<UInt> & Mesh::getUIntData(const ElementType & el_type,
					      const std::string & data_name,
					      const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  const UIntDataMap & map = uint_data(el_type, ghost_type);
  UIntDataMap::const_iterator it = map.find(data_name);

  if(it == map.end())
    AKANTU_EXCEPTION("No data named " << data_name << " in the mesh " << id);

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline UIntDataMap & Mesh::getUIntDataMap(const ElementType & el_type,
					  const GhostType & ghost_type) {
  // AKANTU_DEBUG_ASSERT(uint_data.exists(el_type, ghost_type),
  // 		      "No UIntDataMap for the type (" << ghost_type << ":" << el_type << ")");
  return uint_data(el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const ElementType & type,
			       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  try {
    const Vector<UInt> & conn = connectivities(type, ghost_type);
    AKANTU_DEBUG_OUT();
    return conn.getSize();
  } catch (...) {
    AKANTU_DEBUG_OUT();
    return 0;
  }
}

/* -------------------------------------------------------------------------- */
inline void Mesh::getBarycenter(UInt element, const ElementType & type,
				Real * barycenter,
				GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt * conn_val = getConnectivity(type, ghost_type).values;
  UInt nb_nodes_per_element = getNbNodesPerElement(type);

  Real local_coord[spatial_dimension * nb_nodes_per_element];

  UInt offset = element * nb_nodes_per_element;
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    memcpy(local_coord + n * spatial_dimension,
	   nodes->storage() + conn_val[offset + n] * spatial_dimension,
	   spatial_dimension*sizeof(Real));
  }

  Math::barycenter(local_coord, nb_nodes_per_element, spatial_dimension, barycenter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodesPerElement(const ElementType & type) {
  UInt nb_nodes_per_element = 0;
#define GET_NB_NODES_PER_ELEMENT(type)					\
  nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_NODES_PER_ELEMENT);
#undef GET_NB_NODES_PER_ELEMENT
  return nb_nodes_per_element;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getP1ElementType(const ElementType & type) {
  ElementType p1_type = _not_defined;
#define GET_P1_TYPE(type)				\
  p1_type = ElementClass<type>::getP1ElementType()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_P1_TYPE);
#undef GET_P1_TYPE
  return p1_type;
}

/* -------------------------------------------------------------------------- */
inline ElementKind Mesh::getKind(const ElementType & type) {
  ElementKind kind = _ek_not_defined;
#define GET_KIND(type)				\
  kind = ElementClass<type>::getKind()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_KIND);
#undef GET_KIND
  return kind;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getSpatialDimension(const ElementType & type) {
  UInt spatial_dimension = 0;
#define GET_SPATIAL_DIMENSION(type)				\
  spatial_dimension = ElementClass<type>::getSpatialDimension()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SPATIAL_DIMENSION);
#undef GET_SPATIAL_DIMENSION

  return spatial_dimension;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getFacetType(const ElementType & type) {
  ElementType surface_type = _not_defined;
#define GET_FACET_TYPE(type)					\
  surface_type = ElementClass<type>::getFacetType()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_FACET_TYPE);
#undef GET_FACET_TYPE

  return surface_type;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbFacetsPerElement(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt n_facet = 0;
#define GET_NB_FACET(type)				\
  n_facet = ElementClass<type>::getNbFacetsPerElement()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_FACET);
#undef GET_NB_FACET

  AKANTU_DEBUG_OUT();
  return n_facet;
}

/* -------------------------------------------------------------------------- */
inline types::Matrix<UInt> Mesh::getFacetLocalConnectivity(const ElementType & type) {
  AKANTU_DEBUG_IN();

  types::Matrix<UInt> mat;

#define GET_FACET_CON(type)						\
  mat = ElementClass<type>::getFacetLocalConnectivityPerElement()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_FACET_CON);
#undef GET_FACET_CON

  AKANTU_DEBUG_OUT();
  return mat;
}

/* -------------------------------------------------------------------------- */
inline types::Matrix<UInt> Mesh::getFacetConnectivity(UInt element,
						      const ElementType & type,
						      const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  types::Matrix<UInt> local_facets = getFacetLocalConnectivity(type);
  types::Matrix<UInt> facets(local_facets.rows(), local_facets.cols());

  const Vector<UInt> & conn = connectivities(type, ghost_type);

  for (UInt f = 0; f < facets.rows(); ++f) {
    for (UInt n = 0; n < facets.cols(); ++n) {
      facets(f, n) = conn(element, local_facets(f, n));
    }
  }

  AKANTU_DEBUG_OUT();
  return facets;
}


/* -------------------------------------------------------------------------- */
template<typename T>
inline void Mesh::extractNodalValuesFromElement(const Vector<T> & nodal_values,
						T * local_coord,
						UInt * connectivity,
						UInt n_nodes,
						UInt nb_degree_of_freedom) const {
  for (UInt n = 0; n < n_nodes; ++n) {
    memcpy(local_coord + n * nb_degree_of_freedom,
	   nodal_values.storage() + connectivity[n] * nb_degree_of_freedom,
	   nb_degree_of_freedom * sizeof(T));
  }
}

/* -------------------------------------------------------------------------- */
#define DECLARE_GET_BOUND(Var, var)                               \
  inline void Mesh::get##Var##Bounds(Real * var) const {          \
    for (UInt i = 0; i < spatial_dimension; ++i) {                \
      var[i] = var##_bounds[i];                                   \
    }                                                             \
  }                                                               \

DECLARE_GET_BOUND(Lower, lower)
DECLARE_GET_BOUND(Upper, upper)

DECLARE_GET_BOUND(LocalLower, local_lower)
DECLARE_GET_BOUND(LocalUpper, local_upper)

#undef DECLARE_GET_BOUND

/* -------------------------------------------------------------------------- */
inline void Mesh::addConnectivityType(const ElementType & type){
  getConnectivityPointer(type);
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isPureGhostNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) == -3 : false;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isLocalOrMasterNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) == -2 || (*nodes_type)(n) == -1 : true;
}


/* -------------------------------------------------------------------------- */
inline bool Mesh::isLocalNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) == -1 : true;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isMasterNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) == -2 : false;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isSlaveNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) >= 0 : false;
}


/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNodeGlobalId(UInt local_id) const {
  return nodes_global_ids ? (*nodes_global_ids)(local_id) : local_id;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbGlobalNodes() const {
  return nodes_global_ids ? nb_global_nodes : nodes->getSize();
}

/* -------------------------------------------------------------------------- */
inline Int Mesh::getNodeType(UInt local_id) const {
  return nodes_type ? (*nodes_type)(local_id) : -1;
}
