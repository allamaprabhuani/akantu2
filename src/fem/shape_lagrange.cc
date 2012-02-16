/**
 * @file   shape_lagrange.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Feb 10 11:45:51 2011
 *
 * @brief  lagrangian shape functions class
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
#include "aka_memory.hh"
#include "mesh.hh"
#include "shape_lagrange.hh"
#include "element_class.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ShapeLagrange::ShapeLagrange(Mesh & mesh,
			     const ID & id,
			     const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id),
  shapes("shapes", id),
  shapes_derivatives("shapes_derivatives", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange::precomputeShapesOnControlPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * natural_coords = control_points(type, ghost_type).values;
  UInt nb_points        = control_points(type, ghost_type).getSize();

  UInt size_of_shapes    = ElementClass<type>::getShapeSize();

  UInt nb_element = mesh->getConnectivity(type, ghost_type).getSize();;

  Vector<Real> & shapes_tmp = shapes.alloc(nb_element*nb_points,
					   size_of_shapes,
					   type,
					   ghost_type);

  Real * shapes_val    = shapes_tmp.values;
  for (UInt elem = 0; elem < nb_element; ++elem) {
    ElementClass<type>::computeShapes(natural_coords,
				      nb_points,shapes_val);

    shapes_val += size_of_shapes*nb_points;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange::precomputeShapeDerivativesOnControlPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt size_of_shapesd      = ElementClass<type>::getShapeDerivativesSize();
  UInt nb_points            = control_points(type, ghost_type).getSize();
  Real * natural_coords     = control_points(type, ghost_type).values;


  UInt * elem_val = mesh->getConnectivity(type, ghost_type).values;;
  UInt nb_element = mesh->getConnectivity(type, ghost_type).getSize();

  Vector<Real> & shapes_derivatives_tmp = shapes_derivatives.alloc(nb_element*nb_points,
								   size_of_shapesd,
								   type,
								   ghost_type);

  Real * shapesd_val = shapes_derivatives_tmp.values;
  Real local_coord[spatial_dimension * nb_nodes_per_element];

  for (UInt elem = 0; elem < nb_element; ++elem) {
    mesh->extractNodalCoordinatesFromElement(local_coord,
					     elem_val+elem*nb_nodes_per_element,
					     nb_nodes_per_element);

    computeShapeDerivativesOnCPointsByElement<type>(spatial_dimension,
						    local_coord,
						    nb_nodes_per_element,
						    natural_coords,
						    nb_points,
						    shapesd_val);

    shapesd_val += size_of_shapesd*nb_points;
  }

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange::interpolateOnControlPoints(const Vector<Real> &in_u,
						  Vector<Real> &out_uq,
						  UInt nb_degree_of_freedom,
						  GhostType ghost_type,
						  const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(shapes.exists(type, ghost_type),
		      "No shapes for the type "
		      << shapes.printType(type, ghost_type));

  UInt nb_element = mesh->getNbElement(type, ghost_type);
  UInt * conn_val = mesh->getConnectivity(type, ghost_type).values;

  const Vector<Real> & shapes_loc = shapes(type, ghost_type);

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_points = control_points(type, ghost_type).getSize();
  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  Real * shape_val = shapes_loc.storage();
  Real * u_val     = in_u.values;
  Real * uq_val    = out_uq.values;

  UInt offset_uq   = out_uq.getNbComponent() * nb_points;

  Real * shape = shape_val;
  Real * u = new Real[nb_nodes_per_element * nb_degree_of_freedom];

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shape     = shape_val + filter_elem_val[el] * size_of_shapes * nb_points;
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      memcpy(u + n * nb_degree_of_freedom,
	     u_val + conn_val[el_offset + n] * nb_degree_of_freedom,
	     nb_degree_of_freedom * sizeof(Real));
    }

    /// @f$ u_q = \textbf{N} * u @f$
    Math::matrix_matrix(nb_points, nb_degree_of_freedom, nb_nodes_per_element,
			shape, u, uq_val);

    uq_val += offset_uq;
    if(filter_elements == NULL) {
      shape += size_of_shapes*nb_points;
    }
  }

  delete [] u;

#undef INIT_VARIABLES
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange::gradientOnControlPoints(const Vector<Real> &in_u,
					       Vector<Real> &out_nablauq,
					       UInt nb_degree_of_freedom,
					       GhostType ghost_type,
					       const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(shapes_derivatives.exists(type, ghost_type),
		      "No shapes derivatives for the type "
		      << shapes_derivatives.printType(type, ghost_type));

  UInt nb_element = mesh->getNbElement(type, ghost_type);
  UInt * conn_val = mesh->getConnectivity(type, ghost_type).storage();

  const Vector<Real> & shapesd_loc = shapes_derivatives(type, ghost_type);

  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt size_of_shapes_derivatives = ElementClass<type>::getShapeDerivativesSize();
  UInt nb_points         = control_points(type, ghost_type).getSize();
  UInt element_dimension = ElementClass<type>::getSpatialDimension();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  Real * shaped_val  = shapesd_loc.storage();
  Real * u_val       = in_u.values;
  Real * nablauq_val = out_nablauq.storage();

  UInt offset_nablauq = nb_degree_of_freedom * element_dimension;
  UInt offset_shaped  = nb_nodes_per_element * element_dimension;

  Real * shaped  = shaped_val;
  Real * u       = new Real[nb_nodes_per_element * nb_degree_of_freedom];

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shaped    = shaped_val  + filter_elem_val[el] * size_of_shapes_derivatives * nb_points;
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      memcpy(u + n * nb_degree_of_freedom,
	     u_val + conn_val[el_offset + n] * nb_degree_of_freedom,
	     nb_degree_of_freedom * sizeof(Real));
    }

    for (UInt q = 0; q < nb_points; ++q) {
      /// @f$\nabla u^t = u^t * B^t@f$
      Math::matrixt_matrix(nb_degree_of_freedom, element_dimension, nb_nodes_per_element,
			   u,
			   shaped,
			   nablauq_val);

      nablauq_val += offset_nablauq;
      shaped      += offset_shaped;
    }
  }

  delete [] u;

#undef INIT_VARIABLES
  AKANTU_DEBUG_OUT();
}

// /* -------------------------------------------------------------------------- */
// template <ElementType type>
// void ShapeLagrange::setControlPointsByType(Vector<Real> & points, GhostType ghost_type){
//   control_points(type, ghost_type) = &points;
// }

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange::fieldTimesShapes(const Vector<Real> & field,
				     Vector<Real> & field_times_shapes,
				     GhostType ghost_type) const {
  field_times_shapes.copy(field);
  field_times_shapes.extendComponentsInterlaced(ElementClass<type>::getShapeSize(), 1);

  UInt nb_element = field_times_shapes.getSize() * ElementClass<type>::getShapeSize();

  Real * field_times_shapes_val = field_times_shapes.values;
  Real * shapes_val = shapes(type, ghost_type).values;

  /// compute @f$ rho * \varphi_i @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    *field_times_shapes_val++ *= *shapes_val++;
  }
}


/* -------------------------------------------------------------------------- */
void ShapeLagrange::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Shapes Lagrange [" << std::endl;
  ShapeFunctions::printself(stream, indent + 1);
  shapes.printself(stream, indent + 1);
  shapes_derivatives.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}


/* -------------------------------------------------------------------------- */
/* template instantiation */
/* -------------------------------------------------------------------------- */

#define INSTANCIATE_TEMPLATE_CLASS(type)				\
  template								\
  void ShapeLagrange::precomputeShapesOnControlPoints<type>		\
  (GhostType ghost_type);						\
									\
  template								\
  void ShapeLagrange::precomputeShapeDerivativesOnControlPoints<type>	\
  (GhostType ghost_type);						\
									\
  template								\
  void   ShapeLagrange::gradientOnControlPoints<type>			\
  (const Vector<Real> &in_u,						\
   Vector<Real> &out_nablauq,						\
   UInt nb_degree_of_freedom,						\
   GhostType ghost_type,						\
   const Vector<UInt> * filter_elements) const;				\
									\
  template								\
  void ShapeLagrange::interpolateOnControlPoints<type>			\
  (const Vector<Real> &in_u,						\
   Vector<Real> &out_uq,						\
   UInt nb_degree_of_freedom,						\
   GhostType ghost_type,						\
   const Vector<UInt> * filter_elements) const;				\
									\
  template								\
  void ShapeLagrange::fieldTimesShapes<type>				\
  (const Vector<Real> & field,						\
   Vector<Real> & field_times_shapes,					\
   GhostType ghost_type) const;

AKANTU_BOOST_REGULAR_ELEMENT_LIST(INSTANCIATE_TEMPLATE_CLASS)
#undef INSTANCIATE_TEMPLATE_CLASS

__END_AKANTU__
