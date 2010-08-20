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
inline void FEM::integrate(const Vector<Real> & f,
			   Real * intf,
			   UInt nb_degre_of_freedom,
			   const ElementType & type,
			   const UInt elem,
			   bool local) const {

  Vector<Real> * jac_loc;
#ifdef AKANTU_USE_MPI
  if(local) {
#endif //AKANTU_USE_MPI
    jac_loc     = shapes_derivatives[type];
#ifdef AKANTU_USE_MPI
  } else {
    jac_loc     = ghost_shapes_derivatives[type];
  }
#endif //AKANTU_USE_MPI

  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  UInt size_of_jacobians    = FEM::getJacobianSize(type);

  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degre_of_freedom * nb_quadrature_points,
		      "The vector f do not have the good number of component.");

  Real * f_val    = f.values + elem * f.getNbComponent();
  Real * jac_val  = jac_loc->values + elem * size_of_jacobians;

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
inline UInt FEM::getShapeSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_size;
#define GET_SHAPE_SIZE(type)				\
  shape_size = ElementClass<type>::getShapeSize()

  switch(type) {
  case _line_1       : { GET_SHAPE_SIZE(_line_1      ); break; }
  case _line_2       : { GET_SHAPE_SIZE(_line_2      ); break; }
  case _triangle_1   : { GET_SHAPE_SIZE(_triangle_1  ); break; }
  case _triangle_2   : { GET_SHAPE_SIZE(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_SHAPE_SIZE(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_SHAPE_SIZE(_tetrahedra_2); break; }
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

  UInt shape_derivatives_size;
#define GET_SHAPE_DERIVATIVES_SIZE(type)				\
  shape_derivatives_size = ElementClass<type>::getShapeDerivativesSize()

  switch(type) {
  case _line_1       : { GET_SHAPE_DERIVATIVES_SIZE(_line_1      ); break; }
  case _line_2       : { GET_SHAPE_DERIVATIVES_SIZE(_line_2      ); break; }
  case _triangle_1   : { GET_SHAPE_DERIVATIVES_SIZE(_triangle_1  ); break; }
  case _triangle_2   : { GET_SHAPE_DERIVATIVES_SIZE(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_SHAPE_DERIVATIVES_SIZE(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_SHAPE_DERIVATIVES_SIZE(_tetrahedra_2); break; }
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
inline UInt FEM::getJacobianSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt jacobian_size;
#define GET_JACOBIAN_SIZE(type)				\
  jacobian_size = ElementClass<type>::getJacobianSize()

  switch(type) {
  case _line_1       : { GET_JACOBIAN_SIZE(_line_1      ); break; }
  case _line_2       : { GET_JACOBIAN_SIZE(_line_2      ); break; }
  case _triangle_1   : { GET_JACOBIAN_SIZE(_triangle_1  ); break; }
  case _triangle_2   : { GET_JACOBIAN_SIZE(_triangle_2  ); break; }
  case _tetrahedra_1 : { GET_JACOBIAN_SIZE(_tetrahedra_1); break; }
  case _tetrahedra_2 : { GET_JACOBIAN_SIZE(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

#undef GET_JACOBIAN_SIZE

  AKANTU_DEBUG_OUT();
  return jacobian_size;
}


/* -------------------------------------------------------------------------- */
inline Real FEM::getElementInradius(Real * coord, const ElementType & type) {
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

