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
inline Vector<Real> & Mesh::getNodes() const {
  return *nodes;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodes() const {
  return nodes->getSize();
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> & Mesh::getConnectivity(ElementType type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(connectivities[type] != NULL,
		      "The mesh " << id << " as no element of kind : "<< type);

  AKANTU_DEBUG_OUT();
  return *connectivities[type];
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getConnectivityPointer(ElementType type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return connectivities[type];
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const ElementType & type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(connectivities[type] != NULL,
		      "The mesh " << id << " as no element of kind : "<< type);

  AKANTU_DEBUG_OUT();
  return connectivities[type]->getSize();
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodesPerElement(const ElementType & type) const {
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
inline UInt Mesh::getNbNodesPerElementP1(const ElementType & type) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element_p1;
#define GET_NB_NODES_PER_ELEMENT_P1(type)				\
  nb_nodes_per_element_p1 = ElementClass<type>::getNbNodesPerElementP1()

  switch(type) {
  case _line_1       : { GET_NB_NODES_PER_ELEMENT_P1(_line_1      ); break; }
  case _line_2       : { GET_NB_NODES_PER_ELEMENT_P1(_line_2      ); break; }
  case _triangle_1   : { GET_NB_NODES_PER_ELEMENT_P1(_triangle_1  ); break; }
  case _triangle_2   : { GET_NB_NODES_PER_ELEMENT_P1(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_NB_NODES_PER_ELEMENT_P1(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_NB_NODES_PER_ELEMENT_P1(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_NB_NODES_PER_ELEMENT_P1

  AKANTU_DEBUG_OUT();
  return nb_nodes_per_element_p1;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getSpatialDimension(const ElementType & type) const {
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
inline const ElementType Mesh::getSurfaceElementType(const ElementType & type) const {
  AKANTU_DEBUG_IN();

  ElementType surface_type;
#define GET_SURFACE_TYPE(type)					\
  surface_type = ElementClass<type>::getSurfaceElementType()

  switch(type) {
  case _line_1       : { GET_SURFACE_TYPE(_line_1      ); break; }
  case _line_2       : { GET_SURFACE_TYPE(_line_2      ); break; }
  case _triangle_1   : { GET_SURFACE_TYPE(_triangle_1  ); break; }
  case _triangle_2   : { GET_SURFACE_TYPE(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_SURFACE_TYPE(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_SURFACE_TYPE(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_SURFACE_TYPE

  AKANTU_DEBUG_OUT();
  return surface_type;
}
