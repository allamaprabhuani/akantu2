/**
 * @file   fem_inline.impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jul 19 12:21:36 2010
 *
 * @brief  Implementation of the inline functions of the FEM Class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
inline void FEM::integrate(const Vector<Real> & f,
			   Real * intf,
			   UInt nb_degre_of_freedom,
			   const ElementType & type,
			   const UInt elem,
			   GhostType ghost_type) const {

  Vector<Real> * jac_loc;

  if(ghost_type == _not_ghost) {
    jac_loc     = shapes_derivatives[type];
  } else {
    jac_loc     = ghost_shapes_derivatives[type];
  }

  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degre_of_freedom ,
		      "The vector f do not have the good number of component.");

  Real * f_val    = f.values + elem * f.getNbComponent();
  Real * jac_val  = jac_loc->values + elem * nb_quadrature_points;

  integrate(f_val, jac_val, intf, nb_degre_of_freedom, nb_quadrature_points);
}


/* -------------------------------------------------------------------------- */
inline void FEM::integrate(Real *f, Real *jac, Real * inte,
			   UInt nb_degre_of_freedom,
			   UInt nb_quadrature_points) const {
  memset(inte, 0, nb_degre_of_freedom * sizeof(Real));

  Real *cjac = jac;
  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    for (UInt dof = 0; dof < nb_degre_of_freedom; ++dof) {
      inte[dof] += *f * *cjac;
      ++f;
    }
    ++cjac;
  }
}

/* -------------------------------------------------------------------------- */
inline Mesh & FEM::getMesh() const {
  return *mesh;
}

/* -------------------------------------------------------------------------- */
inline UInt FEM::getNbQuadraturePoints(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = 0;
#define GET_NB_QUAD_POINTS(type)					\
  nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints()

  switch(type) {
  case _segment_2       : { GET_NB_QUAD_POINTS(_segment_2      ); break; }
  case _segment_3       : { GET_NB_QUAD_POINTS(_segment_3      ); break; }
  case _triangle_3   : { GET_NB_QUAD_POINTS(_triangle_3  ); break; }
  case _triangle_6   : { GET_NB_QUAD_POINTS(_triangle_6  ); break; }
  case _tetrahedron_4 : { GET_NB_QUAD_POINTS(_tetrahedron_4); break; }
  case _tetrahedron_10 : { GET_NB_QUAD_POINTS(_tetrahedron_10); break; }
  case _point:
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
inline UInt FEM::getShapeSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_size = 0;
#define GET_SHAPE_SIZE(type)				\
  shape_size = ElementClass<type>::getShapeSize()

  switch(type) {
  case _segment_2       : { GET_SHAPE_SIZE(_segment_2      ); break; }
  case _segment_3       : { GET_SHAPE_SIZE(_segment_3      ); break; }
  case _triangle_3   : { GET_SHAPE_SIZE(_triangle_3  ); break; }
  case _triangle_6   : { GET_SHAPE_SIZE(_triangle_6  ); break; }
  case _tetrahedron_4 : { GET_SHAPE_SIZE(_tetrahedron_4); break; }
  case _tetrahedron_10 : { GET_SHAPE_SIZE(_tetrahedron_10); break; }
  case _point:
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_SHAPE_SIZE

  AKANTU_DEBUG_OUT();
  return shape_size;
}

/* -------------------------------------------------------------------------- */
inline UInt FEM::getShapeDerivativesSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_derivatives_size = 0;
#define GET_SHAPE_DERIVATIVES_SIZE(type)				\
  shape_derivatives_size = ElementClass<type>::getShapeDerivativesSize()

  switch(type) {
  case _segment_2       : { GET_SHAPE_DERIVATIVES_SIZE(_segment_2      ); break; }
  case _segment_3       : { GET_SHAPE_DERIVATIVES_SIZE(_segment_3      ); break; }
  case _triangle_3   : { GET_SHAPE_DERIVATIVES_SIZE(_triangle_3  ); break; }
  case _triangle_6   : { GET_SHAPE_DERIVATIVES_SIZE(_triangle_6  ); break; }
  case _tetrahedron_4 : { GET_SHAPE_DERIVATIVES_SIZE(_tetrahedron_4); break; }
  case _tetrahedron_10 : { GET_SHAPE_DERIVATIVES_SIZE(_tetrahedron_10); break; }
  case _point:
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_SHAPE_DERIVATIVES_SIZE

  AKANTU_DEBUG_OUT();
  return shape_derivatives_size;
}

/* -------------------------------------------------------------------------- */
inline Real FEM::getElementInradius(Real * coord, const ElementType & type) {
  AKANTU_DEBUG_IN();

  Real inradius = 0;

#define GET_INRADIUS(type)						\
  inradius = ElementClass<type>::getInradius(coord);			\

  switch(type) {
  case _segment_2       : { GET_INRADIUS(_segment_2      ); break; }
  case _segment_3       : { GET_INRADIUS(_segment_3      ); break; }
  case _triangle_3   : { GET_INRADIUS(_triangle_3  ); break; }
  case _triangle_6   : { GET_INRADIUS(_triangle_6  ); break; }
  case _tetrahedron_4 : { GET_INRADIUS(_tetrahedron_4); break; }
  case _tetrahedron_10 : { GET_INRADIUS(_tetrahedron_10); break; }
  case _point:
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
