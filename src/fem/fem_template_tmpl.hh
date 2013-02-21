/**
 * @file   fem_template_tmpl.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Nov  5 17:08:50 2012
 *
 * @brief  Template implementation of FEMTemplate
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

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
FEMTemplate<I, S, kind>::FEMTemplate(Mesh & mesh, UInt spatial_dimension,
				      ID id, MemoryID memory_id) :
  FEM(mesh,spatial_dimension,id,memory_id),
  integrator(mesh, id, memory_id),
  shape_functions(mesh, id, memory_id) { }

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
FEMTemplate<I, S, kind>::~FEMTemplate() { }

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::gradientOnQuadraturePoints(const Vector<Real> &u,
							  Vector<Real> &nablauq,
							  const UInt nb_degree_of_freedom,
							  const ElementType & type,
							  const GhostType & ghost_type,
							  const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  UInt nb_element = mesh->getNbElement(type, ghost_type);
  if(filter_elements != NULL) nb_element = filter_elements->getSize();

  UInt element_dimension = mesh->getSpatialDimension(type);
  UInt nb_points         = shape_functions.getControlPoints(type, ghost_type).cols();

  AKANTU_DEBUG_ASSERT(u.getSize() == mesh->getNbNodes(),
		      "The vector u(" << u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom ,
		      "The vector u(" << u.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(nablauq.getNbComponent()
		      == nb_degree_of_freedom * element_dimension,
		      "The vector nablauq(" << nablauq.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(nablauq.getSize() == nb_element * nb_points,
		      "The vector nablauq(" << nablauq.getID()
		      << ") has not the good size.");
#endif

#define COMPUTE_GRADIENT(type)						\
    if (element_dimension == ElementClass<type>::getSpatialDimension()) \
      shape_functions.template gradientOnControlPoints<type>(u,		\
							     nablauq,	\
							     nb_degree_of_freedom, \
							     ghost_type, \
							     filter_elements);

  if(kind == _ek_regular)
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(COMPUTE_GRADIENT);
  else if(kind == _ek_cohesive)
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_GRADIENT);
  else
    AKANTU_DEBUG_TO_IMPLEMENT();
#undef COMPUTE_GRADIENT

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::initShapeFunctions(const GhostType & ghost_type) {
  initShapeFunctions(mesh->getNodes(), ghost_type);
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::initShapeFunctions(const Vector<Real> & nodes,
						 const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it  = mesh->firstType(element_dimension, ghost_type, kind);
  Mesh::type_iterator end = mesh->lastType(element_dimension, ghost_type, kind);
  for(; it != end; ++it) {
    ElementType type = *it;
    integrator.initIntegrator(nodes, type, ghost_type);
    const types::Matrix<Real> control_points =
      getQuadraturePoints(type, ghost_type);
    shape_functions.initShapeFunctions(nodes, control_points, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::integrate(const Vector<Real> & f,
				       Vector<Real> &intf,
				       UInt nb_degree_of_freedom,
				       const ElementType & type,
				       const GhostType & ghost_type,
				       const Vector<UInt> * filter_elements) const{

#ifndef AKANTU_NDEBUG
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
  integrator.template integrate<type>(f,		    \
				      intf,		    \
				      nb_degree_of_freedom, \
				      ghost_type,	    \
				      filter_elements);

  if(kind == _ek_regular)
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_cohesive)
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_structural)
    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INTEGRATE);
  else
    AKANTU_DEBUG_TO_IMPLEMENT();
#undef INTEGRATE
}

/* -------------------------------------------------------------------------- */


template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
Real FEMTemplate<I, S, kind>::integrate(const Vector<Real> & f,
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

  UInt nb_quadrature_points  = getNbQuadraturePoints(type, ghost_type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
		      "The vector f(" << f.getID()
		      << ") has not the good size. (" << f.getSize() << "!=" << nb_quadrature_points * nb_element << ")");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
		      "The vector f(" << f.getID()
		      << ") has not the good number of component.");
#endif

  Real integral = 0.;

#define INTEGRATE(type)							\
  integral = integrator.template integrate<type>(f,			\
						 ghost_type,		\
						 filter_elements);

  if(kind == _ek_regular)
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_cohesive)
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_structural)
    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INTEGRATE);
  else
    AKANTU_DEBUG_TO_IMPLEMENT();
#undef INTEGRATE

  AKANTU_DEBUG_OUT();
  return integral;
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
Real FEMTemplate<I, S, kind>::integrate(const types::RVector & f,
					 const ElementType & type,
					 UInt index,
					 const GhostType & ghost_type) const{


  Real res = 0.;

#define INTEGRATE(type)							\
  res = integrator.template integrate<type>(f,				\
					    index,			\
					    ghost_type);

  if(kind == _ek_regular)
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_cohesive)
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_structural)
    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INTEGRATE);
  else
    AKANTU_DEBUG_TO_IMPLEMENT();
#undef INTEGRATE

  return res;
}


/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::integrateOnQuadraturePoints(const Vector<Real> & f,
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
  AKANTU_DEBUG_ASSERT(intf.getSize() == nb_element * nb_quadrature_points,
		      "The vector intf(" << intf.getID()
		      << ") has not the good size.");
#endif

#define INTEGRATE(type)							\
  integrator.template integrateOnQuadraturePoints<type>(f,		\
							intf,		\
							nb_degree_of_freedom, \
							ghost_type,	\
							filter_elements);

  if(kind == _ek_regular)
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_cohesive)
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
  else if(kind == _ek_structural)
    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INTEGRATE);
  else
    AKANTU_DEBUG_TO_IMPLEMENT();
#undef INTEGRATE
}


/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::interpolateOnQuadraturePoints(const Vector<Real> &u,
							   Vector<Real> &uq,
							   UInt nb_degree_of_freedom,
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

  UInt nb_points         = shape_functions.getControlPoints(type, ghost_type).cols();

  AKANTU_DEBUG_ASSERT(u.getSize() == mesh->getNbNodes(),
		      "The vector u(" << u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom ,
		      "The vector u(" << u.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(uq.getNbComponent() == nb_degree_of_freedom,
		      "The vector uq(" << uq.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(uq.getSize() == nb_element * nb_points,
		      "The vector uq(" << uq.getID()
		      << ") has not the good size.");
#endif


#define INTERPOLATE(type)						\
  shape_functions.template interpolateOnControlPoints<type>(u,		\
							    uq,		\
							    nb_degree_of_freedom, \
							    ghost_type, \
							    filter_elements);

  if(kind == _ek_regular)
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INTERPOLATE);
  else if(kind == _ek_cohesive)
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTERPOLATE);
  else
    AKANTU_DEBUG_TO_IMPLEMENT();
#undef INTERPOLATE

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::computeNormalsOnControlPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  computeNormalsOnControlPoints(mesh->getNodes(),
				ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::computeNormalsOnControlPoints(const Vector<Real> & field,
							     const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  if (ghost_type == _ghost) { AKANTU_DEBUG_TO_IMPLEMENT(); }

  //  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();

  //allocate the normal arrays
  mesh->initByElementTypeVector(normals_on_quad_points, spatial_dimension, element_dimension);

  //loop over the type to build the normals
  Mesh::type_iterator it  = mesh->firstType(element_dimension, ghost_type, kind);
  Mesh::type_iterator end = mesh->lastType(element_dimension, ghost_type, kind);
  for(; it != end; ++it) {
    Vector<Real> & normals_on_quad = normals_on_quad_points(*it, ghost_type);
    computeNormalsOnControlPoints(field, normals_on_quad, *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::computeNormalsOnControlPoints(const Vector<Real> & field,
							     Vector<Real> & normal,
							     const ElementType & type,
							     const GhostType & ghost_type) const {
#define COMPUTE_NORMALS_ON_QUAD(type)					\
  computeNormalsOnControlPoints<type>(field, normal, ghost_type);

  if(kind == _ek_regular)
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(COMPUTE_NORMALS_ON_QUAD);
  else if(kind == _ek_cohesive)
    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMALS_ON_QUAD);
  else
    AKANTU_DEBUG_TO_IMPLEMENT();

#undef COMPUTE_NORMALS_ON_QUAD
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template<ElementType type>
void FEMTemplate<I, S, kind>::computeNormalsOnControlPoints(const Vector<Real> & field,
							     Vector<Real> & normal,
							     const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  if (ghost_type == _ghost) { AKANTU_DEBUG_TO_IMPLEMENT(); }

  UInt spatial_dimension = mesh->getSpatialDimension();
  UInt nb_nodes_per_element  = Mesh::getNbNodesPerElement(type);
  UInt nb_points = getNbQuadraturePoints(type, ghost_type);

  UInt nb_element = mesh->getConnectivity(type, ghost_type).getSize();
  normal.resize(nb_element * nb_points);
  Vector<Real>::iterator< types::Matrix<Real> > normals_on_quad = normal.begin_reinterpret(spatial_dimension,
											   nb_points,
											   nb_element);
  Vector<Real> f_el(0, spatial_dimension * nb_nodes_per_element);
  FEM::extractNodalToElementField(*mesh, field, f_el, type, ghost_type);

  const types::Matrix<Real> quads =
    integrator. template getQuadraturePoints<type>(ghost_type);

  Vector<Real>::iterator< types::Matrix<Real> > f_it = f_el.begin(spatial_dimension, nb_nodes_per_element);

  for (UInt elem = 0; elem < nb_element; ++elem) {
    ElementClass<type>::computeNormalsOnNaturalCoordinates(quads,
							   *f_it,
							   *normals_on_quad);
    ++normals_on_quad;
    ++f_it;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* Matrix lumping functions                                                   */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::assembleFieldLumped(const Vector<Real> & field_1,
						   UInt nb_degree_of_freedom,
						   Vector<Real> & lumped,
						   const Vector<Int> & equation_number,
						   ElementType type,
						   const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

#define ASSEMBLE_LUMPED(type)					\
  assembleLumpedTemplate<type>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type)

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(ASSEMBLE_LUMPED);;

#undef ASSEMBLE_LUMPED
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEMTemplate<I, S, kind>::assembleFieldMatrix(const Vector<Real> & field_1,
						   UInt nb_degree_of_freedom,
						   SparseMatrix & matrix,
						   ElementType type,
						   const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

#define ASSEMBLE_MATRIX(type)					\
  assembleFieldMatrix<type>(field_1, nb_degree_of_freedom,	\
			    matrix,				\
			    ghost_type)

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(ASSEMBLE_MATRIX);;

#undef ASSEMBLE_MATRIX

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEMTemplate<I, S, kind>::assembleLumpedTemplate(const Vector<Real> & field_1,
						      UInt nb_degree_of_freedom,
						      Vector<Real> & lumped,
						      const Vector<Int> & equation_number,
						      const GhostType & ghost_type) const {
  this->template assembleLumpedRowSum<type>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
 */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEMTemplate<I, S, kind>::assembleLumpedRowSum(const Vector<Real> & field_1,
						    UInt nb_degree_of_freedom,
						    Vector<Real> & lumped,
						    const Vector<Int> & equation_number,
						    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt shapes_size = ElementClass<type>::getShapeSize();

  Vector<Real> * field_times_shapes = new Vector<Real>(0, 1);//shapes_size);
  shape_functions.template fieldTimesShapes<type>(field_1, *field_times_shapes, ghost_type);

  UInt nb_element = mesh->getNbElement(type, ghost_type);
  Vector<Real> * int_field_times_shapes = new Vector<Real>(nb_element, shapes_size,
							   "inte_rho_x_shapes");
  integrator.template integrate<type>(*field_times_shapes, *int_field_times_shapes,
				      shapes_size, ghost_type, NULL);
  delete field_times_shapes;

  int_field_times_shapes->extendComponentsInterlaced(nb_degree_of_freedom,1);
  assembleVector(*int_field_times_shapes, lumped, equation_number,nb_degree_of_freedom, type, ghost_type);
  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
 */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEMTemplate<I, S, kind>::assembleLumpedDiagonalScaling(const Vector<Real> & field_1,
							     UInt nb_degree_of_freedom,
							     Vector<Real> & lumped,
							     const Vector<Int> & equation_number,
							     const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element_p1 = ElementClass<type>::getNbNodesPerElement();
  UInt nb_nodes_per_element    = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points    = integrator.template getQuadraturePoints<type>(ghost_type).cols();

  UInt nb_element = field_1.getSize() / nb_quadrature_points;

  Real corner_factor = 0;
  Real mid_factor    = 0;

  if(type == _triangle_6) {
    corner_factor = 1./12.;
    mid_factor    = 1./4.;
  }

  if (type == _tetrahedron_10) {
    corner_factor = 1./32.;
    mid_factor    = 7./48.;
  }

  if (type == _quadrangle_8) {
    corner_factor = 1./36.;
    mid_factor    = 8./36.;
  }

  if (nb_element == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  /// compute @f$ \int \rho dV = \rho V @f$ for each element
  Vector<Real> * int_field_1 = new Vector<Real>(field_1.getSize(), 1,
					      "inte_rho_x_1");
  integrator.template integrate<type>(field_1, *int_field_1, 1, ghost_type, NULL);

  /// distribute the mass of the element to the nodes
  Vector<Real> * lumped_per_node = new Vector<Real>(nb_element, nb_nodes_per_element, "mass_per_node");
  Real * int_field_1_val = int_field_1->values;
  Real * lumped_per_node_val = lumped_per_node->values;

  for (UInt e = 0; e < nb_element; ++e) {
    Real lmass = *int_field_1_val * corner_factor;
    for (UInt n = 0; n < nb_nodes_per_element_p1; ++n)
      *lumped_per_node_val++ = lmass; /// corner points

    lmass = *int_field_1_val * mid_factor;
    for (UInt n = nb_nodes_per_element_p1; n < nb_nodes_per_element; ++n)
      *lumped_per_node_val++ = lmass; /// mid points

    int_field_1_val++;
  }
  delete int_field_1;

  lumped_per_node->extendComponentsInterlaced(nb_degree_of_freedom,1);
  assembleVector(*lumped_per_node, lumped, equation_number, nb_degree_of_freedom, type, ghost_type);
  delete lumped_per_node;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
 */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEMTemplate<I, S, kind>::assembleFieldMatrix(const Vector<Real> & field_1,
						   UInt nb_degree_of_freedom,
						   SparseMatrix & matrix,
						   const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt vect_size   = field_1.getSize();
  UInt shapes_size = ElementClass<type>::getShapeSize();
  UInt lmat_size   = nb_degree_of_freedom * shapes_size;

  const Vector<Real> & shapes = shape_functions.getShapes(type,ghost_type);
  Vector<Real> * modified_shapes = new Vector<Real>(vect_size, lmat_size * nb_degree_of_freedom);
  modified_shapes->clear();
  Vector<Real> * local_mat = new Vector<Real>(vect_size, lmat_size * lmat_size);

  Vector<Real>::iterator<types::RMatrix> shape_vect  = modified_shapes->begin(nb_degree_of_freedom, lmat_size);
  Real * sh  = shapes.values;
  for(UInt q = 0; q < vect_size; ++q) {
    Real * msh = shape_vect->storage();
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      Real * msh_tmp = msh + d * (lmat_size + 1);
      for (UInt s = 0; s < shapes_size; ++s) {
	*msh_tmp = sh[s];
	msh_tmp += nb_degree_of_freedom;
      }
    }
    ++shape_vect;
    sh += shapes_size;
  }

  shape_vect  = modified_shapes->begin(nb_degree_of_freedom, lmat_size);
  Vector<Real>::iterator<types::RMatrix> lmat = local_mat->begin(lmat_size, lmat_size);
  Real * field_val = field_1.values;

  for(UInt q = 0; q < vect_size; ++q) {
    (*lmat).mul<true, false>(*shape_vect, *shape_vect, *field_val);
    ++lmat; ++shape_vect; ++field_val;
  }

  delete modified_shapes;

  UInt nb_element = mesh->getNbElement(type, ghost_type);
  Vector<Real> * int_field_times_shapes = new Vector<Real>(nb_element, lmat_size * lmat_size,
							   "inte_rho_x_shapes");
  integrator.template integrate<type>(*local_mat, *int_field_times_shapes,
				      lmat_size * lmat_size, ghost_type, NULL);
  delete local_mat;

  assembleMatrix(*int_field_times_shapes, matrix, nb_degree_of_freedom, type, ghost_type);
  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline void FEMTemplate<I, S, kind>::inverseMap(const types::Vector<Real> & real_coords,
						UInt element,
						const ElementType & type,
						types::Vector<Real> & natural_coords,
						const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

#define INVERSE_MAP(type) \
  shape_functions.template inverseMap<type>(real_coords, element, natural_coords, ghost_type);

  AKANTU_BOOST_ELEMENT_SWITCH(INVERSE_MAP, AKANTU_NOT_STRUCTURAL_ELEMENT_TYPE);
#undef INVERSE_MAP

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline bool FEMTemplate<I, S, kind>::contains(const types::RVector & real_coords,
					       UInt element,
					       const ElementType & type,
					       const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

  bool contain = false;

#define CONTAINS(type)							\
  contain = shape_functions.template contains<type>(real_coords, element, ghost_type);

  AKANTU_BOOST_ELEMENT_SWITCH(CONTAINS, AKANTU_NOT_STRUCTURAL_ELEMENT_TYPE);

#undef CONTAINS

  AKANTU_DEBUG_OUT();
  return contain;
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline void FEMTemplate<I, S, kind>::computeShapes(const types::RVector & real_coords,
						    UInt element,
						    const ElementType & type,
						    types::RVector & shapes,
						    const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

#define COMPUTE_SHAPES(type) \
  shape_functions.template computeShapes<type>(real_coords,element,shapes,ghost_type);

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(COMPUTE_SHAPES);

#undef COMPUTE_SHAPES

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline UInt FEMTemplate<I, S, kind>::getNbQuadraturePoints(const ElementType & type,
							    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quad_points = 0;

#define GET_NB_QUAD(type)						\
  nb_quad_points =							\
    integrator. template getQuadraturePoints<type>(ghost_type).cols();

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_QUAD);
#undef GET_NB_QUAD

  AKANTU_DEBUG_OUT();
  return nb_quad_points;
}


/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline const Vector<Real> & FEMTemplate<I, S, kind>::getShapes(const ElementType & type,
								const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapes(type, ghost_type));

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline const Vector<Real> & FEMTemplate<I, S, kind>::getShapesDerivatives(const ElementType & type,
									   const GhostType & ghost_type,
									   __attribute__((unused)) UInt id) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapesDerivatives(type, ghost_type));

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline const types::Matrix<Real> &
FEMTemplate<I, S, kind>::getQuadraturePoints(const ElementType & type,
					     const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  const types::Matrix<Real> * ret = NULL;

#define GET_QUADS(type)						\
  ret = &(integrator. template getQuadraturePoints<type>(ghost_type));

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_QUADS);
#undef GET_QUADS

    AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
__END_AKANTU__
#include "shape_lagrange.hh"
#include "integrator_gauss.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEMTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
assembleLumpedTemplate<_triangle_6>(const Vector<Real> & field_1,
				    UInt nb_degree_of_freedom,
				    Vector<Real> & lumped,
				    const Vector<Int> & equation_number,
				    const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_triangle_6>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEMTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
assembleLumpedTemplate<_tetrahedron_10>(const Vector<Real> & field_1,
					UInt nb_degree_of_freedom,
					Vector<Real> & lumped,
					const Vector<Int> & equation_number,
					const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_tetrahedron_10>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void
FEMTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
assembleLumpedTemplate<_quadrangle_8>(const Vector<Real> & field_1,
				      UInt nb_degree_of_freedom,
				      Vector<Real> & lumped,
				      const Vector<Int> & equation_number,
				      const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_quadrangle_8>(field_1, nb_degree_of_freedom,lumped, equation_number, ghost_type);
}


/* -------------------------------------------------------------------------- */
/* Shape Linked specialization                                                */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
__END_AKANTU__
#include "shape_linked.hh"
__BEGIN_AKANTU__

#if defined(AKANTU_STRUCTURAL_MECHANICS)
/* -------------------------------------------------------------------------- */
template <>
inline const Vector<Real> &
FEMTemplate<IntegratorGauss, ShapeLinked, _ek_structural>::getShapesDerivatives(const ElementType & type,
										const GhostType & ghost_type,
										UInt id) const {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapesDerivatives(type, ghost_type, id));

  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}
#endif
