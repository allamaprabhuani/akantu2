/**
 * @file   fem_template.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Fri Feb 11 11:37:47 2011
 *
 * @brief  implementation of the generic FEMTemplate class
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
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "fem.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

template <typename Integ, typename Shape>
FEMTemplate<Integ,Shape>::FEMTemplate(Mesh & mesh, UInt spatial_dimension,
				    FEMID id,MemoryID memory_id)
  :FEM(mesh,spatial_dimension,id,memory_id),
   integrator(mesh),
   shape_functions(mesh)
{
}
/* -------------------------------------------------------------------------- */


template <typename Integ, typename Shape>
FEMTemplate<Integ,Shape>::~FEMTemplate()
{
}

/* -------------------------------------------------------------------------- */

template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::gradientOnQuadraturePoints(const Vector<Real> &u,
							  Vector<Real> &nablauq,
							  const UInt nb_degree_of_freedom,
							  const ElementType & type,
							  GhostType ghost_type,
							  const Vector<UInt> * filter_elements){
  AKANTU_DEBUG_IN();
  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin();
      it != type_list.end();
      ++it) {
    
    ElementType type = *it;
    
#define COMPUTE_GRADIENT(type)						\
    if (element_dimension == ElementClass<type>::getSpatialDimension()) \
      shape_functions.template gradientOnControlPoints<type>(u,		\
							     nablauq,	\
							     nb_degree_of_freedom, \
							     ghost_type, \
							     filter_elements);
    AKANTU_BOOST_ELEMENT_SWITCH(COMPUTE_GRADIENT)
#undef COMPUTE_GRADIENT
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::initShapeFunctions(GhostType ghost_type){
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh->getSpatialDimension();
  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  
  for(it = type_list.begin();
      it != type_list.end();
      ++it) {
    
    ElementType type = *it;
    
#define INIT_SHAPE_FUNCTIONS(type)					\
    if (element_dimension != ElementClass<type>::getSpatialDimension()) continue; \
    integrator.template computeQuadraturePoints<type>();		\
    integrator.template precomputeJacobiansOnQuadraturePoints<type>(spatial_dimension,ghost_type); \
    Vector<Real> & control_points = integrator.template getQuadraturePoints<type>(); \
    shape_functions.template setControlPointsByType<type>(control_points); \
    shape_functions.template precomputeShapesOnControlPoints<type>(ghost_type); \
    if (element_dimension == spatial_dimension)					\
      shape_functions.template precomputeShapeDerivativesOnControlPoints<type>(ghost_type); 
    
    AKANTU_BOOST_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS)
#undef INIT_SHAPE_FUNCTIONS
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::integrate(const Vector<Real> & f,
				       Vector<Real> &intf,
				       UInt nb_degree_of_freedom,
				       const ElementType & type,
				       GhostType ghost_type,
				       const Vector<UInt> * filter_elements) const{

#define INTEGRATE(type)			 \
  integrator.template integrate<type>(f,    \
				      intf,		    \
				      nb_degree_of_freedom, \
				      ghost_type,	    \
				      filter_elements); 
  
    AKANTU_BOOST_ELEMENT_SWITCH(INTEGRATE)
#undef INTEGRATE
}

/* -------------------------------------------------------------------------- */


template <typename Integ, typename Shape>
Real FEMTemplate<Integ,Shape>::integrate(const Vector<Real> & f,
				       const ElementType & type,
				       GhostType ghost_type,
				       const Vector<UInt> * filter_elements) const{
  AKANTU_DEBUG_IN();
#define INTEGRATE(type)						\
  return integrator.template integrate<type>(f,			\
					     ghost_type,	\
					     filter_elements); 
  
  AKANTU_BOOST_ELEMENT_SWITCH(INTEGRATE)
#undef INTEGRATE
    
  return 0;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::interpolateOnQuadraturePoints(const Vector<Real> &u,
							   Vector<Real> &uq,
							   UInt nb_degree_of_freedom,
							   const ElementType & type,
							   GhostType ghost_type,
							   const Vector<UInt> * filter_elements) const{

  AKANTU_DEBUG_IN();
#define INTERPOLATE(type)						\
  return shape_functions.template interpolateOnControlPoints<type>(u,	\
								   uq,	\
								   nb_degree_of_freedom, \
								   ghost_type, \
								   filter_elements); 
  
  AKANTU_BOOST_ELEMENT_SWITCH(INTERPOLATE)
#undef INTERPOLATE
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::computeNormalsOnControlPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();

  //allocate the normal arrays
  if (ghost_type == _not_ghost)
    mesh->initByElementTypeRealVector(normals_on_quad_points,spatial_dimension,element_dimension,
				id,"normals_onquad",ghost_type);
  else{
    AKANTU_DEBUG_ERROR("to be implemented");
  }

  //loop over the type to build the normals
  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin();
      it != type_list.end();
      ++it) {

    ElementType type = *it;

    UInt element_type_spatial_dimension = Mesh::getSpatialDimension(type);

    if(element_type_spatial_dimension != element_dimension) continue;

    UInt nb_nodes_per_element           = Mesh::getNbNodesPerElement(type);
    UInt nb_quad_points = getQuadraturePoints(type).getSize();
    UInt * elem_val;
    UInt nb_element;

    Real * normals_on_quad_val    = NULL;

    if(ghost_type == _not_ghost) {
      elem_val   = mesh->getConnectivity(type).values;
      nb_element = mesh->getConnectivity(type).getSize();
      normals_on_quad_points[type]->resize(nb_element * nb_quad_points);
      normals_on_quad_val =  normals_on_quad_points[type]->values;
    } else {
      elem_val   = mesh->getGhostConnectivity(type).values;
      nb_element = mesh->getGhostConnectivity(type).getSize();
    }

    

    /* ---------------------------------------------------------------------- */
#define COMPUTE_NORMALS_ON_QUAD(type)					\
    do {								\
      Vector<Real> & quads = integrator. template getQuadraturePoints<type>(); \
      UInt nb_points = quads.getSize();					\
      Real local_coord[spatial_dimension * nb_nodes_per_element];	\
      for (UInt elem = 0; elem < nb_element; ++elem) {			\
	mesh->extractNodalCoordinatesFromElement(local_coord,		\
        coord,elem_val+elem*nb_nodes_per_element,nb_nodes_per_element);	\
	ElementClass<type>::computeNormalsOnQuadPoint(local_coord,      \
						  spatial_dimension,	\
						  normals_on_quad_val);	\
	normals_on_quad_val += spatial_dimension*nb_points;      	\
      }									\
    } while(0)
    /* ---------------------------------------------------------------------- */
    
    AKANTU_BOOST_ELEMENT_SWITCH(COMPUTE_NORMALS_ON_QUAD)
#undef COMPUTE_NORMALS_ON_QUAD

  }
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* compatibility functions */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline UInt FEMTemplate<Integ,Shape>::getNbQuadraturePoints(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt nb_quad_points = 0;

#define GET_NB_QUAD(type)						\
  nb_quad_points = integrator. template getQuadraturePoints<type>().getSize();

  AKANTU_BOOST_ELEMENT_SWITCH(GET_NB_QUAD)
#undef GET_NB_QUAD
    
  AKANTU_DEBUG_OUT();
  return nb_quad_points;
}


/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getShapes(const ElementType & type) {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapes(type));

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPES)
#undef GET_SHAPES
    
  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getGhostShapes(const ElementType & type) {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getGhostShapes(type));

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPES)
#undef GET_SHAPES
    
  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getShapesDerivatives(const ElementType & type) {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapesDerivatives(type));

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPES)
#undef GET_SHAPES
    
  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getGhostShapesDerivatives(const ElementType & type) {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getGhostShapesDerivatives(type));

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPES)
#undef GET_SHAPES
    
  AKANTU_DEBUG_OUT();
  return *ret;
}


/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getQuadraturePoints(const ElementType & type) {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;
  
#define GET_QUADS(type)						\
  ret = &(integrator. template getQuadraturePoints<type>());
  
  AKANTU_BOOST_ELEMENT_SWITCH(GET_QUADS)
#undef GET_QUADS
    
    AKANTU_DEBUG_OUT();
  return *ret;
}



/* -------------------------------------------------------------------------- */
/* template instanciation */
/* -------------------------------------------------------------------------- */

template class FEMTemplate<IntegratorGauss,ShapeLagrange>;

__END_AKANTU__
