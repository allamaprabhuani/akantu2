/**
 * @file   fem_template_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date   Tue May 22 20:45:51 2012
 *
 * @brief  FEMTemplate inline implementation
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
template <typename Integ, typename Shape>
inline void FEMTemplate<Integ,Shape>::inverseMap(const types::RVector & real_coords,
						 UInt element,
						 const ElementType & type,
						 types::RVector & natural_coords,
						 const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

#define INVERSE_MAP(type) \
  shape_functions.template inverseMap<type>(real_coords,element,natural_coords,ghost_type);

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INVERSE_MAP);

#undef INVERSE_MAP

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline bool FEMTemplate<Integ,Shape>::contains(const types::RVector & real_coords,
					       UInt element,
					       const ElementType & type,
					       const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

  bool contain = false;

#define CONTAINS(type)							\
  contain = shape_functions.template contains<type>(real_coords,element,ghost_type);

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(CONTAINS);

#undef CONTAINS

  AKANTU_DEBUG_OUT();
  return contain;
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline void FEMTemplate<Integ,Shape>::computeShapes(const types::RVector & real_coords,
						    UInt element,
						    const ElementType & type,
						    types::RVector & shapes,
						    const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

#define COMPUTE_SHAPES(type) \
  shape_functions.template computeShapes<type>(real_coords,element,shapes,ghost_type);

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(COMPUTE_SHAPES);

#undef COMPUTE_SHAPES

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline UInt FEMTemplate<Integ,Shape>::getNbQuadraturePoints(const ElementType & type,
							    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quad_points = 0;

#define GET_NB_QUAD(type)						\
  nb_quad_points =							\
    integrator. template getQuadraturePoints<type>(ghost_type).getSize();

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_NB_QUAD);
#undef GET_NB_QUAD

  AKANTU_DEBUG_OUT();
  return nb_quad_points;
}


/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getShapes(const ElementType & type,
								const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapes(type, ghost_type));

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getShapesDerivatives(const ElementType & type,
									   const GhostType & ghost_type,
									   __attribute__((unused)) UInt id) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapesDerivatives(type, ghost_type));

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getQuadraturePoints(const ElementType & type,
									  const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_QUADS(type)						\
  ret = &(integrator. template getQuadraturePoints<type>(ghost_type));

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_QUADS);
#undef GET_QUADS

    AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
/* Shape Linked specialization                                                */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
template <>
inline bool FEMTemplate<IntegratorGauss,ShapeLinked >
::contains(__attribute__((unused)) const types::RVector & real_coords,
	   __attribute__((unused)) UInt element,
	   __attribute__((unused)) const ElementType & type,
	   __attribute__((unused)) const GhostType & ghost_type) const{

  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline void FEMTemplate<IntegratorGauss,ShapeLinked >
::computeShapes(__attribute__((unused)) const types::RVector & real_coords,
		__attribute__((unused)) UInt element,
		__attribute__((unused)) const ElementType & type,
		__attribute__((unused)) types::RVector & shapes,
		__attribute__((unused)) const GhostType & ghost_type) const{
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline void FEMTemplate<IntegratorGauss,ShapeLinked >
::inverseMap(__attribute__((unused)) const types::RVector & real_coords,
	     __attribute__((unused)) UInt element,
	     __attribute__((unused)) const ElementType & type,
	     __attribute__((unused)) types::RVector & natural_coords,
	     __attribute__((unused)) const GhostType & ghost_type) const{

  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline const Vector<Real> &
FEMTemplate<IntegratorGauss, ShapeLinked>::getShapesDerivatives(const ElementType & type,
								const GhostType & ghost_type,
								UInt id) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapesDerivatives(type, ghost_type, id));

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}
