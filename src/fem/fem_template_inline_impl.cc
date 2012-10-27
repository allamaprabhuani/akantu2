/**
 * @file   fem_template_inline_impl.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Tue May 22 11:21:40 2012
 *
 * @brief  
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

template <>
inline void FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >
::inverseMap(__attribute__((unused)) const types::RVector & real_coords,
	     __attribute__((unused)) UInt element,
	     __attribute__((unused)) const ElementType & type,
	     __attribute__((unused)) types::RVector & natural_coords,
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

template <>
inline bool FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >
::contains(__attribute__((unused)) const types::RVector & real_coords,
	   __attribute__((unused)) UInt element,
	   __attribute__((unused)) const ElementType & type,
	   __attribute__((unused)) const GhostType & ghost_type) const{
  
  AKANTU_DEBUG_TO_IMPLEMENT();
}

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

template <>
inline void FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >
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
::computeShapes(__attribute__((unused)) const types::RVector & real_coords,
		__attribute__((unused)) UInt element,
		__attribute__((unused)) const ElementType & type,
		__attribute__((unused)) types::RVector & shapes,
		__attribute__((unused)) const GhostType & ghost_type) const{
  
  AKANTU_DEBUG_TO_IMPLEMENT();
}

