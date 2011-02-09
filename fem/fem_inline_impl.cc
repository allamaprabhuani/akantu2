/**
 * @file   fem_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Mon Jul 19 12:21:36 2010
 *
 * @brief  Implementation of the inline functions of the FEM Class
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
inline void FEM::integrate(const Vector<Real> & f,
			   Real * intf,
			   UInt nb_degre_of_freedom,
			   const Element & elem,
			   GhostType ghost_type) const {

  Vector<Real> * jac_loc;

  if(ghost_type == _not_ghost) {
    jac_loc     = jacobians[elem.type];
  } else {
    jac_loc     = ghost_jacobians[elem.type];
  }

  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(elem.type);
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degre_of_freedom ,
		      "The vector f do not have the good number of component.");

  Real * f_val    = f.values + elem.element * f.getNbComponent();
  Real * jac_val  = jac_loc->values + elem.element * nb_quadrature_points;

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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_NB_QUAD_POINTS)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPE_SIZE)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_SIZE)
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_INRADIUS)
#undef GET_INRADIUS

  AKANTU_DEBUG_OUT();
  return inradius;
}

/* -------------------------------------------------------------------------- */
