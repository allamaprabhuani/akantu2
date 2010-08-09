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
    volume = ElementClass<type>::getVolume(local_coord);		\
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
inline void FEM::integrate(const Vector<Real> & f,
			   Real * intf,
			   UInt nb_degre_of_freedom,
			   const ElementType & type,
			   const UInt elem) {

  UInt nb_quadrature_points = getNbQuadraturePoints(type);
  UInt size_of_jacobians    = jacobians[type]->getNbComponent();

  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degre_of_freedom * nb_quadrature_points,
		      "The vector f do not have the good number of component.");

  Real * f_val    = f.values + elem * f.getNbComponent();
  Real * jac_val  = jacobians[type]->values + elem * size_of_jacobians;

  integrate(f_val, jac_val, intf, nb_degre_of_freedom, nb_quadrature_points);
}


/* -------------------------------------------------------------------------- */
inline void FEM::integrate(Real *f, Real *jac, Real * inte,
			   UInt nb_degre_of_freedom,
			   UInt nb_quadrature_points){
  //  for (UInt n = 0; n < nb_nodes_per_element; ++n) {

  memset(inte, 0, nb_degre_of_freedom * sizeof(Real));

  Real *cjac = jac;
  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    for (UInt dof = 0; dof < nb_degre_of_freedom; ++dof) {
      inte[dof] += *f * *cjac;
      ++f;
    }
    ++cjac;
  }
    //    inte += nb_degre_of_freedom;
    //  }
}

/* -------------------------------------------------------------------------- */
inline Mesh & FEM::getMesh() const {
  return *mesh;
}

/* -------------------------------------------------------------------------- */
inline const Mesh::ConnectivityTypeList & FEM::getConnectivityTypeList() const {
  return mesh->getTypeList();
}

/* -------------------------------------------------------------------------- */
inline UInt FEM::getNbNodes() const {
  return mesh->getNbNodes();
}

/* -------------------------------------------------------------------------- */
inline UInt FEM::getNbElement(const ElementType & type) const {
  return mesh->getConnectivity(type).getSize();
}

/* -------------------------------------------------------------------------- */
inline UInt FEM::getSpatialDimension(const ElementType & type) const {
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
inline UInt FEM::getNbQuadraturePoints(const ElementType & type) const {
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

/* -------------------------------------------------------------------------- */
inline UInt FEM::getNbNodesPerElement(const ElementType & type) const {
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
inline Real FEM::getElementInradius(Real * coord, const ElementType & type) const {
  AKANTU_DEBUG_IN();

  Real inradius;

#define GET_INRADIUS(type)						\
  inradius = ElementClass<type>::getInradius(coord);			\

  switch(type) {
  case _line_1       : { GET_INRADIUS(_line_1      ); break; }
  case _line_2       : { GET_INRADIUS(_line_2      ); break; }
  case _triangle_1   : { GET_INRADIUS(_triangle_1  ); break; }
  case _triangle_2   : { GET_INRADIUS(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_INRADIUS(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_INRADIUS(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_INRADIUS

  AKANTU_DEBUG_OUT();
  return inradius;
}


/* -------------------------------------------------------------------------- */
inline const Vector<Real> & FEM::getShapes(const ElementType & type) const {
  AKANTU_DEBUG_ASSERT(shapes[type] != NULL,
		      "No shapes of the type : " << type << " in " << id);
  return *shapes[type];
}

/* -------------------------------------------------------------------------- */
inline const Vector<Real> & FEM::getShapesDerivatives(const ElementType & type) const {
  AKANTU_DEBUG_ASSERT(shapes_derivatives[type] != NULL,
		      "No shapes derivatives of the type : " << type << " in " << id);
  return *shapes_derivatives[type];
}

/* -------------------------------------------------------------------------- */

