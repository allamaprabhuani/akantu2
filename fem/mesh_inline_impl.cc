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

#ifdef AKANTU_USE_MPI
/* -------------------------------------------------------------------------- */
inline UInt Mesh::ghostElementToLinearized(const Element & elem) {
  AKANTU_DEBUG_ASSERT(elem.type < _max_element_type &&
		      elem.element < ghost_types_offsets.values[elem.type+1],
		      "The ghost element " << elem
		      << "does not exists in the mesh " << id);

  return ghost_types_offsets.values[elem.types] + elem.element;
}

/* -------------------------------------------------------------------------- */
inline Element Mesh::ghostLinearizedToElement (UInt linearized_element) {
  UInt t;
  for (t = _not_defined + 1;
       linearized_element > ghost_types_offsets.values[t] && t <= _max_element_type; ++t);

  AKANTU_DEBUG_ASSERT(t < _max_element_type,
		      "The ghost linearized element " << linearized_element
		      << "does not exists in the mesh " << id);

  t--;
  return Element((ElementType) t, linearized_element - ghost_types_offsets[t]);
}

/* -------------------------------------------------------------------------- */
inline void Mesh::updateGhostTypesOffsets() {
  UInt count = 0;
  for (ElementType t = _not_defined;  t <= _max_element_type; ++t) {
    ghost_types_offsets.values[t] = count;
    count += (t == _max_element_type || ghost_connectivities[t] == NULL) ?
      0 : ghost_connectivities[t]->getSize();
  }
}
#endif //AKANTU_USE_MPI

/* -------------------------------------------------------------------------- */
inline const Mesh::ConnectivityTypeList & Mesh::getConnectivityTypeList(bool local) const {
  #ifdef AKANTU_USE_MPI
    if(local) {
#endif //AKANTU_USE_MPI
      return type_set;
#ifdef AKANTU_USE_MPI
    } else {
      return ghost_type_set;
    }
#endif //AKANTU_USE_MPI
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getConnectivityPointer(ElementType type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return connectivities[type];
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
		      "The mesh " << id << " has no element of kind : "<< type);

  AKANTU_DEBUG_OUT();
  return connectivities[type]->getSize();
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
inline ElementType Mesh::getElementP1(const ElementType & type) {
  AKANTU_DEBUG_IN();

  ElementType element_p1;
#define GET_ELEMENT_P1(type)				\
  element_p1 = ElementClass<type>::getElementP1()

  switch(type) {
  case _line_1       : { GET_ELEMENT_P1(_line_1      ); break; }
  case _line_2       : { GET_ELEMENT_P1(_line_2      ); break; }
  case _triangle_1   : { GET_ELEMENT_P1(_triangle_1  ); break; }
  case _triangle_2   : { GET_ELEMENT_P1(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_ELEMENT_P1(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_ELEMENT_P1(_tetrahedra_2); break; }
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
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_SURFACE_TYPE

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
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_SURFACE_TYPE

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
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_SURFACE_TYPE
  
  AKANTU_DEBUG_OUT();
  return facet_conn;
}

