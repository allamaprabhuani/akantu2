/**
 * @file   mesh_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 14 23:58:08 2010
 *
 * @brief  Implementation of the inline functions of the mesh class
 *
 * @section LICENSE
 *
 * <insert license here>
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
inline Vector<Real> * Mesh::getNormalsPointer(ElementType type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return normals[type];
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
inline UInt Mesh::getNbNodesPerElement(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element;
#define GET_NB_NODES_PER_ELEMENT(type)					\
  nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement()

  switch(type) {
  case _line_1       : { GET_NB_NODES_PER_ELEMENT(_line_1      ); break; }
  case _line_2       : { GET_NB_NODES_PER_ELEMENT(_line_2      ); break; }
  case _triangle_1   : { GET_NB_NODES_PER_ELEMENT(_triangle_1  ); break; }
  case _triangle_2   : { GET_NB_NODES_PER_ELEMENT(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_NB_NODES_PER_ELEMENT(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_NB_NODES_PER_ELEMENT(_tetrahedra_2); break; }
  case _point        : { GET_NB_NODES_PER_ELEMENT(_point       ); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_NB_NODES_PER_ELEMENT

  AKANTU_DEBUG_OUT();
  return nb_nodes_per_element;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getP1ElementType(const ElementType & type) {
  AKANTU_DEBUG_IN();

  ElementType element_p1;
#define GET_ELEMENT_P1(type)				\
  element_p1 = ElementClass<type>::getP1ElementType()

  switch(type) {
  case _line_1       : { GET_ELEMENT_P1(_line_1      ); break; }
  case _line_2       : { GET_ELEMENT_P1(_line_2      ); break; }
  case _triangle_1   : { GET_ELEMENT_P1(_triangle_1  ); break; }
  case _triangle_2   : { GET_ELEMENT_P1(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_ELEMENT_P1(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_ELEMENT_P1(_tetrahedra_2); break; }
  case _point        : { GET_ELEMENT_P1(_point       ); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_NB_NODES_PER_ELEMENT_P1

  AKANTU_DEBUG_OUT();
  return element_p1;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getSpatialDimension(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension;
#define GET_SPATIAL_DIMENSION(type)					\
  spatial_dimension = ElementClass<type>::getSpatialDimension()

  switch(type) {
  case _line_1       : { GET_SPATIAL_DIMENSION(_line_1      ); break; }
  case _line_2       : { GET_SPATIAL_DIMENSION(_line_2      ); break; }
  case _triangle_1   : { GET_SPATIAL_DIMENSION(_triangle_1  ); break; }
  case _triangle_2   : { GET_SPATIAL_DIMENSION(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_SPATIAL_DIMENSION(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_SPATIAL_DIMENSION(_tetrahedra_2); break; }
  case _point        : { GET_SPATIAL_DIMENSION(_point       ); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_SPATIAL_DIMENSION

  AKANTU_DEBUG_OUT();
  return spatial_dimension;
}

/* -------------------------------------------------------------------------- */
inline const ElementType Mesh::getFacetElementType(const ElementType & type) {
  AKANTU_DEBUG_IN();

  ElementType surface_type;
#define GET_FACET_TYPE(type)					\
  surface_type = ElementClass<type>::getFacetElementType()

  switch(type) {
  case _line_1       : { GET_FACET_TYPE(_line_1      ); break; }
  case _line_2       : { GET_FACET_TYPE(_line_2      ); break; }
  case _triangle_1   : { GET_FACET_TYPE(_triangle_1  ); break; }
  case _triangle_2   : { GET_FACET_TYPE(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_FACET_TYPE(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_FACET_TYPE(_tetrahedra_2); break; }
  case _point        : { GET_FACET_TYPE(_point       ); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_FACET_TYPE

  AKANTU_DEBUG_OUT();
  return surface_type;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbFacetsPerElementType(const ElementType & type) const {
  AKANTU_DEBUG_IN();

  UInt n_facet;
#define GET_NB_FACET(type)					\
  n_facet = ElementClass<type>::getNbFacetsPerElement()

  switch(type) {
  case _line_1       : { GET_NB_FACET(_line_1      ); break; }
  case _line_2       : { GET_NB_FACET(_line_2      ); break; }
  case _triangle_1   : { GET_NB_FACET(_triangle_1  ); break; }
  case _triangle_2   : { GET_NB_FACET(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_NB_FACET(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_NB_FACET(_tetrahedra_2); break; }
  case _point        : { GET_NB_FACET(_point       ); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_NB_FACET

  AKANTU_DEBUG_OUT();
  return n_facet;
}


/* -------------------------------------------------------------------------- */
inline UInt ** Mesh::getFacetLocalConnectivityPerElementType(const ElementType & type) const {
  AKANTU_DEBUG_IN();

  UInt ** facet_conn;
#define GET_FACET_CON(type)                                      \
  facet_conn = ElementClass<type>::getFacetLocalConnectivityPerElement()

  switch(type) {
  case _line_1       : { GET_FACET_CON(_line_1      ); break; }
  case _line_2       : { GET_FACET_CON(_line_2      ); break; }
  case _triangle_1   : { GET_FACET_CON(_triangle_1  ); break; }
  case _triangle_2   : { GET_FACET_CON(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_FACET_CON(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_FACET_CON(_tetrahedra_2); break; }
  case _point        : { GET_FACET_CON(_point       ); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_FACET_CON

  AKANTU_DEBUG_OUT();
  return facet_conn;
}
