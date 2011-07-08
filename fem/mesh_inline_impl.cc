/**
 * @file   mesh_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 14 23:58:08 2010
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

/* -------------------------------------------------------------------------- */
inline UInt Mesh::elementToLinearized(const Element & elem) {
  AKANTU_DEBUG_ASSERT(elem.type < _max_element_type &&
		      elem.element < types_offsets.values[elem.type+1],
		      "The element " << elem
		      << "does not exists in the mesh " << id);

  return types_offsets.values[elem.type] + elem.element;
}

/* -------------------------------------------------------------------------- */
inline Element Mesh::linearizedToElement (UInt linearized_element) {
  UInt t;
  for (t = _not_defined + 1;
       linearized_element > types_offsets.values[t] && t <= _max_element_type; ++t);

  AKANTU_DEBUG_ASSERT(t < _max_element_type,
		      "The linearized element " << linearized_element
		      << "does not exists in the mesh " << id);

  return Element((ElementType) t, linearized_element-types_offsets.values[t]);
}

/* -------------------------------------------------------------------------- */
inline void Mesh::updateTypesOffsets(GhostType ghost_type) {
  types_offsets.clear();
  ConnectivityTypeList::const_iterator it;
  for (it = type_set.begin(); it != type_set.end(); ++it)
    types_offsets(*it) = connectivities(*it, ghost_type)->getSize();

  for (UInt t = _not_defined + 1;  t <= _max_element_type; ++t)
    types_offsets(t) += types_offsets(t - 1);
}

// /* -------------------------------------------------------------------------- */
// inline UInt Mesh::ghostElementToLinearized(const Element & elem) {
//   AKANTU_DEBUG_ASSERT(elem.type < _max_element_type &&
// 		      elem.element < ghost_types_offsets.values[elem.type+1],
// 		      "The ghost element " << elem
// 		      << "does not exists in the mesh " << id);

//   return ghost_types_offsets.values[elem.type] +
//     elem.element +
//     types_offsets.values[_max_element_type];
// }

// /* -------------------------------------------------------------------------- */
// inline Element Mesh::ghostLinearizedToElement (UInt linearized_element) {
//   AKANTU_DEBUG_ASSERT(linearized_element >= types_offsets.values[_max_element_type],
// 		      "The linearized element " << linearized_element
// 		      << "is not a ghost element in the mesh " << id);


//   linearized_element -= types_offsets.values[_max_element_type];
//   UInt t;
//   for (t = _not_defined + 1;
//        linearized_element > ghost_types_offsets.values[t] && t <= _max_element_type; ++t);

//   AKANTU_DEBUG_ASSERT(t < _max_element_type,
// 		      "The ghost linearized element " << linearized_element
// 		      << "does not exists in the mesh " << id);

//   t--;
//   return Element((ElementType) t, linearized_element - ghost_types_offsets.values[t]);
// }

// /* -------------------------------------------------------------------------- */
// inline void Mesh::updateGhostTypesOffsets() {
//   ghost_types_offsets.clear();
//   ByElementTypeUInt::const_iterator it;
//   for (it = ghost_connectivities.begin(); it != ghost_connectivities.end(); ++it)
//     ghost_types_offsets(it->first) = it->second->getSize();

//   for (UInt t = _not_defined + 1;  t <= _max_element_type; ++t)
//     ghost_types_offsets(t) += ghost_types_offsets(t - 1);
// }

/* -------------------------------------------------------------------------- */
inline const Mesh::ConnectivityTypeList & Mesh::getConnectivityTypeList(GhostType ghost_type) const {
  if(ghost_type == _not_ghost)
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

  if(!connectivities.exists(type, ghost_type)) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    std::stringstream sstr;
    sstr << id << ":";
    if (ghost_type == _ghost) sstr << "ghost_";
    sstr << "connectivity:" << type;
    connectivities(type, ghost_type) = &(alloc<UInt>(sstr.str(),
						     0,
						     nb_nodes_per_element));

    AKANTU_DEBUG_INFO("The connectivity vector for the type "
		      << type << " created");

    if (ghost_type == _not_ghost) type_set.insert(type);
    else ghost_type_set.insert(type);

    updateTypesOffsets(ghost_type);
  }

  AKANTU_DEBUG_OUT();
  return connectivities(type, ghost_type);
}

// /* -------------------------------------------------------------------------- */
// inline Vector<UInt> * Mesh::getGhostConnectivityPointer(const ElementType & type) {
//   AKANTU_DEBUG_IN();

//   if(ghost_connectivities[type] == NULL) {
//     UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

//     std::stringstream sstr;
//     sstr << id << ":ghost_connectivity:" << type;
//     ghost_connectivities[type] = &(alloc<UInt>(sstr.str(),
// 					 0,
// 					 nb_nodes_per_element));
//     ghost_type_set.insert(type);

//     AKANTU_DEBUG_INFO("The connectivity vector for the type "
// 		      << type << " created");

//     updateGhostTypesOffsets();
//   }

//   AKANTU_DEBUG_OUT();
//   return ghost_connectivities[type];
// }

/* -------------------------------------------------------------------------- */
inline const Mesh & Mesh::getInternalFacetsMesh() const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  if (!internal_facets_mesh)
    AKANTU_DEBUG_ERROR("Internal facets mesh was not created before access "
		       << "=> use mesh utils to that purpose");
  return *internal_facets_mesh;
}

/* -------------------------------------------------------------------------- */
inline Mesh * Mesh::getInternalFacetsMeshPointer() {
  AKANTU_DEBUG_IN();

  if (!internal_facets_mesh){
    std::stringstream name(this->id);
    name << ":internalfacets";
    internal_facets_mesh = new Mesh(this->spatial_dimension-1,*this->nodes,name.str());
  }

  AKANTU_DEBUG_OUT();
  return internal_facets_mesh;
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getSurfaceIDPointer(const ElementType & type, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  if(!surface_id.exists(type, ghost_type)) {
    std::stringstream sstr;
    sstr << id << ":surface_id:" << type;
    surface_id(type, ghost_type) = &(alloc<UInt>(sstr.str(),
						 0,
						 1));

    AKANTU_DEBUG_INFO("The surface id vector for the type "
		      << type << " created");
  }

  AKANTU_DEBUG_OUT();
  return surface_id(type, ghost_type);
}

/* -------------------------------------------------------------------------- */

// inline Vector<UInt> * Mesh::getReversedElementsPBCPointer(const ElementType & type) {
//   AKANTU_DEBUG_IN();

//   if(reversed_elements_pbc[type] == NULL) {
//     AKANTU_DEBUG_ERROR("There are no reversed elements for the type" << type);
//   }

//   AKANTU_DEBUG_OUT();
//   return reversed_elements_pbc[type];
// }


/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getUIntDataPointer(const ElementType & el_type,
					       const std::string & data_name,
					       const GhostType & ghost_type) {
  //  AKANTU_DEBUG_IN();

  Vector<UInt> * data;
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

  AKANTU_DEBUG_ASSERT(it != map.end(),
		      "No data named " << data_name << " in the mesh " << id);

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline UIntDataMap & Mesh::getUIntDataMap(const ElementType & el_type,
					  const GhostType & ghost_type) {
  // AKANTU_DEBUG_ASSERT(uint_data.exists(el_type, ghost_type),
  // 		      "No UIntDataMap for the type (" << ghost_type << ":" << el_type << ")");
  return uint_data(el_type, ghost_type);
};

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const ElementType & type,
			       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  const ByElementTypeUInt & const_conn = connectivities;

  AKANTU_DEBUG_OUT();
  return const_conn(type, ghost_type)->getSize();
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
	   nodes->values + conn_val[offset + n] * spatial_dimension,
	   spatial_dimension*sizeof(Real));
  }

  Math::barycenter(local_coord, nb_nodes_per_element, spatial_dimension, barycenter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void Mesh::setSurfaceIDsFromIntData(std::string & data_name) {
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const Mesh::ConnectivityTypeList & type_list = getConnectivityTypeList(gt);
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension - 1) continue;

      UIntDataMap & map = uint_data(*it, gt);
      UIntDataMap::iterator it_data = map.find(data_name);
      AKANTU_DEBUG_ASSERT(it_data != map.end(),
			  "No data named " << data_name
			  << " present in the mesh " << id
			  << " for the element type " << *it);
      AKANTU_DEBUG_ASSERT(!surface_id.exists(*it, gt),
			  "Surface id for type (" << gt << ":" << *it
			  << ") already set to the vector " << surface_id(*it, gt)->getID());

      surface_id(*it, gt) = it_data->second;
    }
  }
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodesPerElement(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = 0;
#define GET_NB_NODES_PER_ELEMENT(type)					\
  nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_NB_NODES_PER_ELEMENT);
#undef GET_NB_NODES_PER_ELEMENT

  AKANTU_DEBUG_OUT();
  return nb_nodes_per_element;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getP1ElementType(const ElementType & type) {
  AKANTU_DEBUG_IN();

  ElementType element_p1 = _not_defined;
#define GET_ELEMENT_P1(type)				\
  element_p1 = ElementClass<type>::getP1ElementType()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_ELEMENT_P1);
#undef GET_NB_NODES_PER_ELEMENT_P1

  AKANTU_DEBUG_OUT();
  return element_p1;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getSpatialDimension(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = 0;
#define GET_SPATIAL_DIMENSION(type)					\
  spatial_dimension = ElementClass<type>::getSpatialDimension()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SPATIAL_DIMENSION);
#undef GET_SPATIAL_DIMENSION

  AKANTU_DEBUG_OUT();
  return spatial_dimension;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getFacetElementType(const ElementType & type) {
  AKANTU_DEBUG_IN();

  ElementType surface_type = _not_defined;
#define GET_FACET_TYPE(type)					\
  surface_type = ElementClass<type>::getFacetElementType()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_FACET_TYPE);
#undef GET_FACET_TYPE

  AKANTU_DEBUG_OUT();
  return surface_type;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbFacetsPerElement(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt n_facet = 0;
#define GET_NB_FACET(type)					\
  n_facet = ElementClass<type>::getNbFacetsPerElement()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_NB_FACET);
#undef GET_NB_FACET

  AKANTU_DEBUG_OUT();
  return n_facet;
}


/* -------------------------------------------------------------------------- */
inline UInt ** Mesh::getFacetLocalConnectivity(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt ** facet_conn = NULL;
#define GET_FACET_CON(type)                                      \
  facet_conn = ElementClass<type>::getFacetLocalConnectivityPerElement()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_FACET_CON);
#undef GET_FACET_CON

  AKANTU_DEBUG_OUT();
  return facet_conn;
}

/* -------------------------------------------------------------------------- */
inline void Mesh::extractNodalCoordinatesFromElement(Real * local_coord,
						     UInt * connectivity,
						     UInt n_nodes){
  for (UInt n = 0; n < n_nodes; ++n) {
    memcpy(local_coord + n * spatial_dimension,
	   nodes->values + connectivity[n] * spatial_dimension,
	   spatial_dimension * sizeof(Real));
  }
}
/* -------------------------------------------------------------------------- */
// inline void Mesh::extractNodalCoordinatesFromPBCElement(Real * local_coord,
// 							UInt * connectivity,
// 							UInt n_nodes){

//   // get the min max of the element coordinates in all directions
//   Real min[3];
//   Real max[3];

//   for (UInt k = 0; k < spatial_dimension; k++) {
//     min[k] = std::numeric_limits<double>::max();
//     max[k] = std::numeric_limits<double>::min();
//   }

//   for (UInt nd = 0; nd < n_nodes; nd++) {
//     Real * coord = nodes->values + connectivity[nd] * spatial_dimension;
//     for (UInt k = 0; k < spatial_dimension; ++k) {
//       min[k] = std::min(min[k],coord[k]);
//       max[k] = std::max(max[k],coord[k]);
//     }
//   }
//   Real center[3];
//   // compute the center of the element
//   for (UInt k = 0; k < spatial_dimension; ++k) {
//     center[k] = (max[k] + min[k])/2;
//   }
//   // reverse the coordinates that needs it
//   Real * lcoord = local_coord;
//   for (UInt n = 0; n < n_nodes; ++n) {
//     Real * coord = nodes->values + connectivity[n] * spatial_dimension;
//     for (UInt k = 0; k < spatial_dimension; ++k, ++coord,++lcoord) {
//       // if not at a border then copy normally the node
//       if (!pbc_directions[k] || fabs(xmax[k] - *coord) > Math::tolerance){
// 	*lcoord = *coord;
//       }
//       // else if distance from center is larger than global box size
//       // reverting is needed
//       else if (fabs(min[k] - *coord) > size[k]/2-Math::tolerance){
// 	*lcoord = *coord - size[k];
//       }
//       else {
// 	*lcoord = *coord;
//       }
//     }
//   }
// }

/* -------------------------------------------------------------------------- */
inline void Mesh::addConnecticityType(const ElementType & type){
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

/* -------------------------------------------------------------------------- */
/* ByElementType                                                              */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<class Stored>
inline std::string ByElementType<Stored>::printType(const ElementType & type,
						    const GhostType & ghost_type) {
  std::stringstream sstr; sstr << "(" << ghost_type << ":" << type << ")";
  return sstr.str();
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline bool ByElementType<Stored>::exists(ElementType type, GhostType ghost_type) const {
  return this->getData(ghost_type).find(type) != this->getData(ghost_type).end();
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline const Stored & ByElementType<Stored>::operator()(const ElementType & type,
						const GhostType & ghost_type) const {
  typename ByElementType<Stored>::DataMap::const_iterator it =
    this->getData(ghost_type).find(type);

  AKANTU_DEBUG_ASSERT(it != this->getData(ghost_type).end(),
		      "No element of type "
		      << ByElementType<Stored>::printType(type, ghost_type)
		      << " in this ByElementType<"
		      << debug::demangle(typeid(Stored).name()) << "> class");
  return it->second;
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline Stored & ByElementType<Stored>::operator()(const ElementType & type,
					  const GhostType & ghost_type) {
  typename ByElementType<Stored>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end()) {
    ByElementType<Stored>::DataMap & data = this->getData(ghost_type);
    const std::pair<typename DataMap::iterator, bool> & res =
      data.insert(std::pair<ElementType, Stored>(type, Stored()));
    it = res.first;
  }

  return it->second;
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline typename ByElementType<Stored>::DataMap & ByElementType<Stored>::getData(GhostType ghost_type) {
  if(ghost_type == _not_ghost) return data;
  else if (ghost_type == _ghost) return ghost_data;
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline const typename ByElementType<Stored>::DataMap & ByElementType<Stored>::getData(GhostType ghost_type) const {
  if(ghost_type == _not_ghost) return data;
  else if (ghost_type == _ghost) return ghost_data;
}

/* -------------------------------------------------------------------------- */
/// Works only if stored is a pointer to a class with a printself method
template<class Stored>
void ByElementType<Stored>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << "ByElementType<" << debug::demangle(typeid(Stored).name()) << "> [" << std::endl;
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const ByElementType<Stored>::DataMap & data = getData(gt);
    typename ByElementType<Stored>::DataMap::const_iterator it;
    for(it = data.begin(); it != data.end(); ++it) {
      stream << space << debug::demangle(typeid(Stored).name())
	     << ByElementType<Stored>::printType(it->first, gt) << " [" << std::endl;
      it->second->printself(stream, indent + 2);
    }
  }
}

/* -------------------------------------------------------------------------- */
template<class Stored>
ByElementType<Stored>::ByElementType() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Stored>
ByElementType<Stored>::~ByElementType() {

}
