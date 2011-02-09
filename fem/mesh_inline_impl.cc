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
inline void Mesh::updateTypesOffsets() {
  UInt count = 0;
  for (UInt t = _not_defined;  t <= _max_element_type; ++t) {
    types_offsets.values[t] = count;
    count += (t == _max_element_type || connectivities[t] == NULL) ?
      0 : connectivities[t]->getSize();
  }
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::ghostElementToLinearized(const Element & elem) {
  AKANTU_DEBUG_ASSERT(elem.type < _max_element_type &&
		      elem.element < ghost_types_offsets.values[elem.type+1],
		      "The ghost element " << elem
		      << "does not exists in the mesh " << id);

  return ghost_types_offsets.values[elem.type] +
    elem.element +
    types_offsets.values[_max_element_type];
}

/* -------------------------------------------------------------------------- */
inline Element Mesh::ghostLinearizedToElement (UInt linearized_element) {
  AKANTU_DEBUG_ASSERT(linearized_element >= types_offsets.values[_max_element_type],
		      "The linearized element " << linearized_element
		      << "is not a ghost element in the mesh " << id);


  linearized_element -= types_offsets.values[_max_element_type];
  UInt t;
  for (t = _not_defined + 1;
       linearized_element > ghost_types_offsets.values[t] && t <= _max_element_type; ++t);

  AKANTU_DEBUG_ASSERT(t < _max_element_type,
		      "The ghost linearized element " << linearized_element
		      << "does not exists in the mesh " << id);

  t--;
  return Element((ElementType) t, linearized_element - ghost_types_offsets.values[t]);
}

/* -------------------------------------------------------------------------- */
inline void Mesh::updateGhostTypesOffsets() {
  UInt count = 0;
  for (UInt t = _not_defined;  t <= _max_element_type; ++t) {
    ghost_types_offsets.values[t] = count;
    count += (t == _max_element_type || ghost_connectivities[t] == NULL) ?
      0 : ghost_connectivities[t]->getSize();
  }
}

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
    std::stringstream sstr;
    sstr << id << ":nodes_global_ids";
    nodes_global_ids = &(alloc<UInt>(sstr.str(),
				     nodes->getSize(),
				     1));
  }

  AKANTU_DEBUG_OUT();
  return nodes_global_ids;
}


/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getConnectivityPointer(ElementType type) {
  AKANTU_DEBUG_IN();

  if(connectivities[type] == NULL) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    std::stringstream sstr;
    sstr << id << ":connectivity:" << type;
    connectivities[type] = &(alloc<UInt>(sstr.str(),
					 0,
					 nb_nodes_per_element));
    type_set.insert(type);

    AKANTU_DEBUG_INFO("The connectivity vector for the type "
		      << type << " created");

    updateTypesOffsets();
  }

  AKANTU_DEBUG_OUT();
  return connectivities[type];
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getGhostConnectivityPointer(ElementType type) {
  AKANTU_DEBUG_IN();

  if(ghost_connectivities[type] == NULL) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    std::stringstream sstr;
    sstr << id << ":ghost_connectivity:" << type;
    ghost_connectivities[type] = &(alloc<UInt>(sstr.str(),
					 0,
					 nb_nodes_per_element));
    ghost_type_set.insert(type);

    AKANTU_DEBUG_INFO("The connectivity vector for the type "
		      << type << " created");

    updateGhostTypesOffsets();
  }

  AKANTU_DEBUG_OUT();
  return ghost_connectivities[type];
}

/* -------------------------------------------------------------------------- */
inline const Mesh & Mesh::getInternalFacetsMesh() const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  if (!internal_facets_mesh) AKANTU_DEBUG_ERROR("internal facets mesh was not created before access => use mesh utils to that purpose");
  return *internal_facets_mesh;
}

/* -------------------------------------------------------------------------- */
inline Mesh * Mesh::getInternalFacetsMeshPointer() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  if (!internal_facets_mesh){
    std::stringstream name(this->id);
    name << ":internalfacets";
    internal_facets_mesh = new Mesh(this->spatial_dimension-1,*this->nodes,name.str());
  }

  return internal_facets_mesh;
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getSurfaceIdPointer(ElementType type) {
  AKANTU_DEBUG_IN();

  if(surface_id[type] == NULL) {
    std::stringstream sstr;
    sstr << id << ":surface_id:" << type;
    surface_id[type] = &(alloc<UInt>(sstr.str(),
				     0,
				     1));

    AKANTU_DEBUG_INFO("The surface id vector for the type "
		      << type << " created");
  }

  AKANTU_DEBUG_OUT();
  return surface_id[type];
}


/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const ElementType & type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(connectivities[type] != NULL,
		      "No element of kind : " << type << " in " << id);

  AKANTU_DEBUG_OUT();
  return connectivities[type]->getSize();
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbGhostElement(const ElementType & type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(ghost_connectivities[type] != NULL,
		      "No element of kind : " << type << " in " << id);

  AKANTU_DEBUG_OUT();
  return ghost_connectivities[type]->getSize();
}

/* -------------------------------------------------------------------------- */
inline void Mesh::getBarycenter(UInt element, ElementType type, Real * barycenter,
				GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt * conn_val;
  if (ghost_type == _not_ghost) {
    conn_val = getConnectivity(type).values;
  } else {
    conn_val = getGhostConnectivity(type).values;
  }
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
inline const Vector<UInt> & Mesh::getSurfaceId(const ElementType & type) const{
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(surface_id[type] != NULL,
		      "No element of kind : " << type << " in " << id);

  AKANTU_DEBUG_OUT();
  return *surface_id[type];
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodesPerElement(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = 0;
#define GET_NB_NODES_PER_ELEMENT(type)					\
  nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_NB_NODES_PER_ELEMENT)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_ELEMENT_P1)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SPATIAL_DIMENSION)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_FACET_TYPE)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_NB_FACET)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_FACET_CON)
#undef GET_FACET_CON

  AKANTU_DEBUG_OUT();
  return facet_conn;
}
