/**
 * @file   fem_inline.impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jul 19 12:21:36 2010
 *
 * @brief  Implementation of the inline functions of the FEM Class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline Mesh & FEM::getMesh() const {
  return *mesh;
}

/* -------------------------------------------------------------------------- */
inline Real FEM::volume(ElementType type, Int element) {
  AKANTU_DEBUG_IN();

  Real * coord = mesh->getNodes().values;
  Vector<UInt> & conn = mesh->getConnectivity(type);

  UInt * elem_val  = conn.values;
  Real volume = 0;

#define GET_VOLUME(type)						\
  do {									\
    UInt nb_nodes_per_element =						\
      ElementClass<type>::getNbNodesPerElement();			\
    Real local_coord[spatial_dimension * nb_nodes_per_element];		\
    int offset = element * nb_nodes_per_element;			\
    for (UInt id = 0; id < nb_nodes_per_element; ++id) {		\
      memcpy(local_coord + id * spatial_dimension,			\
	     coord + elem_val[offset + id] * spatial_dimension,		\
	     spatial_dimension*sizeof(Real));				\
    }									\
    volume = ElementClass<type>::volume(local_coord);			\
  } while(0)

  switch(type) {
  case _line_1       : { GET_VOLUME(_line_1      ); break; }
  case _line_2       : { GET_VOLUME(_line_2      ); break; }
  case _triangle_1   : { GET_VOLUME(_triangle_1  ); break; }
  case _triangle_2   : { GET_VOLUME(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_VOLUME(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_VOLUME(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_VOLUME
  AKANTU_DEBUG_OUT();
  return volume;
}

/* -------------------------------------------------------------------------- */
inline UInt FEM::getNbQuadraturePoints(ElementType type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points;
#define GET_NB_QUAD_POINTS(type)					\
  nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints()

  switch(type) {
  case _line_1       : { GET_NB_QUAD_POINTS(_line_1      ); break; }
  case _line_2       : { GET_NB_QUAD_POINTS(_line_2      ); break; }
  case _triangle_1   : { GET_NB_QUAD_POINTS(_triangle_1  ); break; }
  case _triangle_2   : { GET_NB_QUAD_POINTS(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_NB_QUAD_POINTS(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_NB_QUAD_POINTS(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_NB_QUAD_POINTS

  AKANTU_DEBUG_OUT();
  return nb_quadrature_points;
}
