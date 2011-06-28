/**
 * @file   fem_template.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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

#define COMPUTE_GRADIENT(type)						\
    if (element_dimension == ElementClass<type>::getSpatialDimension()) \
      shape_functions.template gradientOnControlPoints<type>(u,		\
							     nablauq,	\
							     nb_degree_of_freedom, \
							     ghost_type, \
							     filter_elements);
    AKANTU_BOOST_ELEMENT_SWITCH(COMPUTE_GRADIENT);
#undef COMPUTE_GRADIENT

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
    if (element_dimension != ElementClass<type>::getSpatialDimension()) \
      continue;								\
    integrator.template computeQuadraturePoints<type>();		\
    integrator.								\
      template precomputeJacobiansOnQuadraturePoints<type>(ghost_type);	\
    integrator.								\
      template checkJacobians<type>(ghost_type);			\
    Vector<Real> & control_points =					\
      integrator.template getQuadraturePoints<type>();			\
    shape_functions.							\
      template setControlPointsByType<type>(control_points);		\
    shape_functions.							\
      template precomputeShapesOnControlPoints<type>(ghost_type);	\
    if (element_dimension == spatial_dimension)				\
      shape_functions.							\
	template precomputeShapeDerivativesOnControlPoints<type>(ghost_type);
    
    AKANTU_BOOST_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
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

    AKANTU_BOOST_ELEMENT_SWITCH(INTEGRATE);
#undef INTEGRATE
}

/* -------------------------------------------------------------------------- */


template <typename Integ, typename Shape>
Real FEMTemplate<Integ,Shape>::integrate(const Vector<Real> & f,
				       const ElementType & type,
				       GhostType ghost_type,
				       const Vector<UInt> * filter_elements) const{
  AKANTU_DEBUG_IN();

  Real integral = 0.;

#define INTEGRATE(type)						\
  integral = integrator.template integrate<type>(f,			\
						 ghost_type,		\
						 filter_elements);

  AKANTU_BOOST_ELEMENT_SWITCH(INTEGRATE);
#undef INTEGRATE

  AKANTU_DEBUG_OUT();
  return integral;
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
  shape_functions.template interpolateOnControlPoints<type>(u,		\
							    uq,		\
							    nb_degree_of_freedom, \
							    ghost_type, \
							    filter_elements);

  AKANTU_BOOST_ELEMENT_SWITCH(INTERPOLATE);
#undef INTERPOLATE

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::computeNormalsOnControlPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh->getNodes().values;
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
    UInt * elem_val = mesh->getConnectivity(type,ghost_type).values;
    UInt nb_element = mesh->getConnectivity(type,ghost_type).getSize();

    Real * normals_on_quad_val    = NULL;

    if(ghost_type == _not_ghost) {
      normals_on_quad_points[type]->resize(nb_element * nb_quad_points);
      normals_on_quad_val =  normals_on_quad_points[type]->values;
    } else {
      //TODO something should be done here ?
    }

    /* ---------------------------------------------------------------------- */
#define COMPUTE_NORMALS_ON_QUAD(type)					\
    do {								\
      Vector<Real> & quads = integrator. template getQuadraturePoints<type>(); \
      UInt nb_points = quads.getSize();					\
      Real local_coord[spatial_dimension * nb_nodes_per_element];	\
      for (UInt elem = 0; elem < nb_element; ++elem) {			\
	mesh->extractNodalCoordinatesFromElement(local_coord,		\
						 elem_val+elem*nb_nodes_per_element, \
						 nb_nodes_per_element);	\
	ElementClass<type>::computeNormalsOnQuadPoint(local_coord,      \
						  spatial_dimension,	\
						  normals_on_quad_val);	\
	normals_on_quad_val += spatial_dimension*nb_points;	 	\
      }									\
    } while(0)
    /* ---------------------------------------------------------------------- */

    AKANTU_BOOST_ELEMENT_SWITCH(COMPUTE_NORMALS_ON_QUAD);
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_NB_QUAD);
#undef GET_NB_QUAD

  AKANTU_DEBUG_OUT();
  return nb_quad_points;
}


/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getShapes(const ElementType & type,
								const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapes(type,ghost_type));

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
inline const Vector<Real> & FEMTemplate<Integ,Shape>::getShapesDerivatives(const ElementType & type,
									   const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  const Vector<Real> * ret = NULL;

#define GET_SHAPES(type)						\
  ret = &(shape_functions.getShapesDerivatives(type,ghost_type));

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPES);
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

  AKANTU_BOOST_ELEMENT_SWITCH(GET_QUADS);
#undef GET_QUADS

    AKANTU_DEBUG_OUT();
  return *ret;
}


/* -------------------------------------------------------------------------- */
/* Matrix lumping functions                                                   */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::assembleFieldLumped(const Vector<Real> & field_1,
						   UInt nb_degree_of_freedom,
						   Vector<Real> & lumped,
						   const Vector<Int> & equation_number,
						   ElementType type,
						   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#define ASSEMBLE_LUMPED(type)					\
  assembleLumpedTemplate<type>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type)

  AKANTU_BOOST_ELEMENT_SWITCH(ASSEMBLE_LUMPED);;

#undef ASSEMBLE_LUMPED
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
void FEMTemplate<Integ,Shape>::assembleFieldMatrix(const Vector<Real> & field_1,
						   UInt nb_degree_of_freedom,
						   SparseMatrix & matrix,
						   ElementType type,
						   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#define ASSEMBLE_MATRIX(type)					\
  assembleFieldMatrix<type>(field_1, nb_degree_of_freedom,	\
			    matrix,				\
			    ghost_type)

  AKANTU_BOOST_ELEMENT_SWITCH(ASSEMBLE_MATRIX);;

#undef ASSEMBLE_MATRIX

  AKANTU_DEBUG_OUT();
}




/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
template <ElementType type>
void FEMTemplate<Integ,Shape>::assembleLumpedTemplate(const Vector<Real> & field_1,
						      UInt nb_degree_of_freedom,
						      Vector<Real> & lumped,
						      const Vector<Int> & equation_number,
						      GhostType ghost_type) {
  this->template assembleLumpedRowSum<type>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
//template <typename Integ, typename Shape>
template <>
template <>
void FEMTemplate<IntegratorGauss,ShapeLagrange>::
assembleLumpedTemplate<_triangle_6>(const Vector<Real> & field_1,
				    UInt nb_degree_of_freedom,
				    Vector<Real> & lumped,
				    const Vector<Int> & equation_number,
				    GhostType ghost_type) {
  assembleLumpedDiagonalScaling<_triangle_6>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
void FEMTemplate<IntegratorGauss,ShapeLagrange>::
assembleLumpedTemplate<_tetrahedron_10>(const Vector<Real> & field_1,
					UInt nb_degree_of_freedom,
					Vector<Real> & lumped,
					const Vector<Int> & equation_number,
					GhostType ghost_type) {
  assembleLumpedDiagonalScaling<_tetrahedron_10>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
void FEMTemplate<IntegratorGauss,ShapeLagrange>::assembleLumpedTemplate<_quadrangle_8>(const Vector<Real> & field_1,
										       UInt nb_degree_of_freedom,
										       Vector<Real> & lumped,
										       const Vector<Int> & equation_number,
										       GhostType ghost_type) {
  assembleLumpedDiagonalScaling<_quadrangle_8>(field_1, nb_degree_of_freedom,lumped, equation_number, ghost_type);
}


/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
 */
template <typename Integ, typename Shape>
template <ElementType type>
void FEMTemplate<Integ,Shape>::assembleLumpedRowSum(const Vector<Real> & field_1,
						    UInt nb_degree_of_freedom,
						    Vector<Real> & lumped,
						    const Vector<Int> & equation_number,
						    
						    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt shapes_size = ElementClass<type>::getShapeSize();

  Vector<Real> * field_times_shapes = new Vector<Real>(0, 1);//shapes_size);
  shape_functions.template fieldTimesShapes<type>(field_1, *field_times_shapes, ghost_type);

  Vector<Real> * int_field_times_shapes = new Vector<Real>(0, shapes_size,
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
template <typename Integ, typename Shape>
template <ElementType type>
void FEMTemplate<Integ,Shape>::assembleLumpedDiagonalScaling(const Vector<Real> & field_1,
							     UInt nb_degree_of_freedom,
							     Vector<Real> & lumped,
							     const Vector<Int> & equation_number,
							     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element_p1 = Mesh::getNbNodesPerElement(Mesh::getP1ElementType(type));
  UInt nb_nodes_per_element    = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points    = integrator.template getQuadraturePoints<type>().getSize();

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
template <typename Integ, typename Shape>
template <ElementType type>
void FEMTemplate<Integ,Shape>::assembleFieldMatrix(const Vector<Real> & field_1,
						   UInt nb_degree_of_freedom,
						   SparseMatrix & matrix,
						   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt vect_size   = field_1.getSize();
  UInt shapes_size = ElementClass<type>::getShapeSize();
  UInt lmat_size   = nb_degree_of_freedom * shapes_size;

  const Vector<Real> & shapes = shape_functions.getShapes(type,ghost_type);
  Vector<Real> * modified_shapes = new Vector<Real>(vect_size, lmat_size * nb_degree_of_freedom);
  modified_shapes->clear();
  Vector<Real> * local_mat = new Vector<Real>(vect_size, lmat_size * lmat_size);

  Vector<Real>::iterator<types::Matrix> shape_vect  = modified_shapes->begin(nb_degree_of_freedom, lmat_size);
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
  Vector<Real>::iterator<types::Matrix> lmat        = local_mat->begin(lmat_size, lmat_size);
  Real * field_val = field_1.values;

  for(UInt q = 0; q < vect_size; ++q) {
    (*lmat).mul<true, false>(*shape_vect, *shape_vect, *field_val);
    ++lmat; ++shape_vect; ++field_val;
  }

  delete modified_shapes;

  Vector<Real> * int_field_times_shapes = new Vector<Real>(0, lmat_size * lmat_size,
							   "inte_rho_x_shapes");
  integrator.template integrate<type>(*local_mat, *int_field_times_shapes,
				      lmat_size * lmat_size, ghost_type, NULL);
  delete local_mat;

  assembleMatrix(*int_field_times_shapes, matrix, nb_degree_of_freedom, type, ghost_type);
  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* template instanciation                                                     */
/* -------------------------------------------------------------------------- */

template class FEMTemplate<IntegratorGauss,ShapeLagrange>;

__END_AKANTU__
