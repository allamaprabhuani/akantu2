/**
 * @file   fem_template_cohesive_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Nov  5 15:05:16 2012
 *
 * @brief  Inline part of FEMTemplate for cohesive elements
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

__END_AKANTU__
#include "shape_cohesive.hh"
#include "integrator_cohesive.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <>
inline void FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >::
inverseMap(__attribute__((unused)) const Vector<Real> & real_coords,
	   __attribute__((unused)) UInt element,
	   __attribute__((unused)) const ElementType & type,
	   __attribute__((unused)) Vector<Real> & natural_coords,
	   __attribute__((unused)) const GhostType & ghost_type) const{
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline bool FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >::
contains(__attribute__((unused)) const Vector<Real> & real_coords,
	 __attribute__((unused)) UInt element,
	 __attribute__((unused)) const ElementType & type,
	 __attribute__((unused)) const GhostType & ghost_type) const{

  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline void FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >::
computeShapes(__attribute__((unused)) const Vector<Real> & real_coords,
	      __attribute__((unused)) UInt element,
	      __attribute__((unused)) const ElementType & type,
	      __attribute__((unused)) Vector<Real> & shapes,
	      __attribute__((unused)) const GhostType & ghost_type) const{
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline UInt FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >::
getNbQuadraturePoints(const ElementType & type,
		      const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quad_points = 0;

#define GET_NB_QUAD(type)						\
  nb_quad_points =							\
    integrator. getQuadraturePoints<type>(ghost_type).cols();

  //  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_NB_QUAD);
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(GET_NB_QUAD);
#undef GET_NB_QUAD

  AKANTU_DEBUG_OUT();
  return nb_quad_points;
}

/* -------------------------------------------------------------------------- */
template <>
inline const Array<Real> & FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >::
getShapes(const ElementType & type,
	  const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const Array<Real> * ret = NULL;

#define GET_SHAPES(type)				\
  ret = &(shape_functions.getShapes(type, ghost_type));

  //  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_SHAPES);
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

