/**
 * @file   fem_template_cohesive.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Oct 31 16:24:42 2012
 *
 * @brief  Specialization for cohesive element
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
#include "fem_template.hh"
#include "shape_cohesive.hh"
#include "integrator_cohesive.hh"

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
inline bool FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >
::contains(__attribute__((unused)) const types::RVector & real_coords,
	   __attribute__((unused)) UInt element,
	   __attribute__((unused)) const ElementType & type,
	   __attribute__((unused)) const GhostType & ghost_type) const{

  AKANTU_DEBUG_TO_IMPLEMENT();
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
/* -------------------------------------------------------------------------- */
template <>
void FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >::
					       initShapeFunctions(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh->getSpatialDimension();
  Mesh::type_iterator it  = mesh->firstType(element_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator end = mesh->lastType(element_dimension, ghost_type, _ek_cohesive);
  for(; it != end; ++it) {
    ElementType type = *it;

#define INIT_SHAPE_FUNCTIONS(type)					\
    integrator.computeQuadraturePoints<type>(ghost_type);		\
    integrator.								\
      precomputeJacobiansOnQuadraturePoints<type>(ghost_type);		\
    integrator.								\
      checkJacobians<type>(ghost_type);					\
    const Vector<Real> & control_points =				\
      integrator.getQuadraturePoints<type>(ghost_type);			\
    shape_functions.							\
      setControlPointsByType<type>(control_points, ghost_type);		\
    shape_functions.							\
      precomputeShapesOnControlPoints<type>(ghost_type);		\
    if (element_dimension == spatial_dimension)				\
      shape_functions.							\
     	precomputeShapeDerivativesOnControlPoints<type>(ghost_type);

    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
#undef INIT_SHAPE_FUNCTIONS
  }
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* compatibility functions */
/* -------------------------------------------------------------------------- */
template <>
inline UInt FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >
::getNbQuadraturePoints(const ElementType & type,
			const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quad_points = 0;

#define GET_NB_QUAD(type)						\
  nb_quad_points =							\
    integrator. getQuadraturePoints<type>(ghost_type).getSize();

  //  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_NB_QUAD);
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(GET_NB_QUAD);
#undef GET_NB_QUAD

  AKANTU_DEBUG_OUT();
  return nb_quad_points;
}

/* -------------------------------------------------------------------------- */
template <>
Real FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >::integrate(const Vector<Real> & f,
											       const ElementType & type,
											       const GhostType & ghost_type,
											       const Vector<UInt> * filter_elements) const{
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
//   std::stringstream sstr; sstr << ghost_type;
//   AKANTU_DEBUG_ASSERT(sstr.str() == nablauq.getTag(),
// 		      "The vector " << nablauq.getID() << " is not taged " << ghost_type);
  UInt nb_element = mesh->getNbElement(type, ghost_type);
  if(filter_elements != NULL) nb_element = filter_elements->getSize();

  UInt nb_quadrature_points  = getNbQuadraturePoints(type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
		      "The vector f(" << f.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
		      "The vector f(" << f.getID()
		      << ") has not the good number of component.");
#endif

  Real integral = 0.;

#define INTEGRATE(type)							\
  integral = integrator. integrate<type>(f,				\
						 ghost_type,		\
						 filter_elements);

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
#undef INTEGRATE

  AKANTU_DEBUG_OUT();
  return integral;
}

/* -------------------------------------------------------------------------- */
template <>
void FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >
::integrate(const Vector<Real> & f,
	    Vector<Real> &intf,
	    UInt nb_degree_of_freedom,
	    const ElementType & type,
	    const GhostType & ghost_type,
	    const Vector<UInt> * filter_elements) const{

#ifndef AKANTU_NDEBUG
//   std::stringstream sstr; sstr << ghost_type;
//   AKANTU_DEBUG_ASSERT(sstr.str() == nablauq.getTag(),
// 		      "The vector " << nablauq.getID() << " is not taged " << ghost_type);
  UInt nb_element = mesh->getNbElement(type, ghost_type);
  if(filter_elements != NULL) nb_element = filter_elements->getSize();

  UInt nb_quadrature_points  = getNbQuadraturePoints(type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
		      "The vector f(" << f.getID() << " size " << f.getSize()
		      << ") has not the good size (" << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom ,
		      "The vector f(" << f.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degree_of_freedom,
		      "The vector intf(" << intf.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getSize() == nb_element,
		      "The vector intf(" << intf.getID()
		      << ") has not the good size.");
#endif

#define INTEGRATE(type)					    \
  integrator. integrate<type>(f,			    \
			      intf,			    \
			      nb_degree_of_freedom,	    \
			      ghost_type,		    \
			      filter_elements);

  //    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INTEGRATE);
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
#undef INTEGRATE
}

/* -------------------------------------------------------------------------- */
template <>
inline const Vector<Real> & FEMTemplate<IntegratorCohesive<IntegratorGauss>,ShapeCohesive<ShapeLagrange> >
::getShapes(const ElementType & type,
	    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapes(type, ghost_type));

  //  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_SHAPES);
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <>
void FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> >::
gradientOnQuadraturePoints(__attribute__((unused)) const Vector<Real> &u,
			   __attribute__((unused)) Vector<Real> &nablauq,
			   __attribute__((unused)) const UInt nb_degree_of_freedom,
			   __attribute__((unused)) const ElementType & type,
			   __attribute__((unused)) const GhostType & ghost_type,
			   __attribute__((unused)) const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template class FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> >;
