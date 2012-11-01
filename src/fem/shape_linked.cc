/**
 * @file   shape_linked.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date   Fri Jul 15 19:41:58 2011
 *
 * @brief  ShapeLinked implementation
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
#include "shape_linked.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ShapeLinked::ShapeLinked(Mesh & mesh, const ID & id, const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id)
{

}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLinked::precomputeShapesOnControlPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();
  UInt nb_nodes_per_element           = Mesh::getNbNodesPerElement(type);

  Real * natural_coords = control_points(type, ghost_type).storage();
  UInt nb_points = control_points(type, ghost_type).getSize();

  UInt size_of_shapes    = ElementClass<type>::getShapeSize();

  UInt * elem_val;
  UInt nb_element;

  std::string ghost = "";
  if(ghost_type == _ghost) {
    ghost = "ghost_";
  }

  elem_val   = mesh->getConnectivity(type, ghost_type).storage();
  nb_element = mesh->getConnectivity(type, ghost_type).getSize();

  UInt nb_shape_functions = ElementClass<type>::getNbShapeFunctions();

  Vector<Real> ** shapes_tmp = new Vector<Real> *[nb_shape_functions];
  for (UInt s = 0; s < nb_shape_functions; ++s) {
    std::stringstream sstr_shapes;
    sstr_shapes << id << ":" << ghost << "shapes:" << type << ":" << s;
    shapes_tmp[s] = &(alloc<Real>(sstr_shapes.str(),
				  nb_element*nb_points,
				  size_of_shapes));

    Real * shapes_val    = shapes_tmp[s]->values;

    Real local_coord[spatial_dimension * nb_nodes_per_element];
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh->extractNodalValuesFromElement(mesh->getNodes(),
					  local_coord,
					  elem_val+elem*nb_nodes_per_element,
					  nb_nodes_per_element,
					  spatial_dimension);

      ElementClass<type>::computeShapes(natural_coords,
					nb_points,
					shapes_val, local_coord, s);

      shapes_val += size_of_shapes*nb_points;
    }
  }

  shapes(type, ghost_type) = shapes_tmp;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLinked::precomputeShapeDerivativesOnControlPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  // Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();

  UInt nb_nodes_per_element           = Mesh::getNbNodesPerElement(type);
  UInt size_of_shapesd   = ElementClass<type>::getShapeDerivativesSize();
  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();
  Real * natural_coords = control_points(type, ghost_type).storage();
  // UInt nb_points = control_points(type, ghost_type)->getSize();

  UInt * elem_val;
  UInt nb_element;
  std::string ghost = "";

  if(ghost_type == _ghost) {
    ghost = "ghost_";
  }

  elem_val   = mesh->getConnectivity(type, ghost_type).values;
  nb_element = mesh->getConnectivity(type, ghost_type).getSize();

  UInt nb_shape_functions = ElementClass<type>::getNbShapeFunctions();

  Vector<Real> ** shapes_derivatives_tmp = new Vector<Real> *[nb_shape_functions];
  for (UInt s = 0; s < nb_shape_functions; ++s) {
    std::stringstream sstr_shapesd;
    sstr_shapesd << id << ":" << ghost << "shapes_derivatives:" << type << ":" << s;
    shapes_derivatives_tmp[s] = &(alloc<Real>(sstr_shapesd.str(),
					      nb_element*nb_quadrature_points,
					      size_of_shapesd));
    Real * shapesd_val   = shapes_derivatives_tmp[s]->values;
    Real local_coord[spatial_dimension * nb_nodes_per_element];
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh->extractNodalValuesFromElement(mesh->getNodes(),
					  local_coord,
					  elem_val+elem*nb_nodes_per_element,
					  nb_nodes_per_element,
					  spatial_dimension);

      // compute shape derivatives
      ElementClass<type>::computeShapeDerivatives(natural_coords,
						  nb_quadrature_points,
						  spatial_dimension,
						  shapesd_val, local_coord, s);

      shapesd_val += size_of_shapesd*nb_quadrature_points;
    }
  }

  shapes_derivatives(type, ghost_type) = shapes_derivatives_tmp;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLinked::interpolateOnControlPoints(const Vector<Real> &in_u,
					     Vector<Real> &out_uq,
					     UInt nb_degree_of_freedom,
					     const GhostType & ghost_type,
					     const Vector<UInt> * filter_elements,
					     bool accumulate,
					     UInt id_shape,
					     UInt num_degre_of_freedom_to_interpolate,
					     UInt num_degre_of_freedom_interpolated) const {
  AKANTU_DEBUG_IN();

  Vector<Real> * shapes_loc;
  UInt nb_element;
  UInt * conn_val;

  shapes_loc = shapes(type, ghost_type)[id_shape];
  nb_element = mesh->getNbElement(type, ghost_type);
  conn_val   = mesh->getConnectivity(type, ghost_type).values;

  AKANTU_DEBUG_ASSERT(shapes_loc != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_points = control_points(type, ghost_type).getSize();
  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  if(!accumulate)
    out_uq.clear();

  Real * shape_val = shapes_loc->values;
  Real * u_val     = in_u.values;
  Real * uq_val    = out_uq.values;

  UInt offset_uq   = out_uq.getNbComponent()*nb_points;

  Real * shape = shape_val;
  Real * u = static_cast<Real *>(calloc(nb_nodes_per_element,
					sizeof(Real)));

  Real * uq = static_cast<Real *>(calloc(nb_points,
					sizeof(Real)));

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shape     = shape_val + filter_elem_val[el] * size_of_shapes*nb_points;
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      u[n] = u_val[conn_val[el_offset + n] * nb_degree_of_freedom + num_degre_of_freedom_to_interpolate];
    }

    /// Uq = Shape * U : matrix product
    Math::matrix_matrix(nb_points, 1, nb_nodes_per_element,
			shape, u, uq);

    for (UInt p = 0; p < nb_points; ++p) {
      uq_val[p*nb_degree_of_freedom + num_degre_of_freedom_interpolated] += uq[p];
    }

   uq_val += offset_uq;
    if(filter_elements == NULL) {
      shape += size_of_shapes*nb_points;
    }
  }

  free(u);
  free(uq);

#undef INIT_VARIABLES
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLinked::gradientOnControlPoints(const Vector<Real> &in_u,
					  Vector<Real> &out_nablauq,
					  UInt nb_degree_of_freedom,
					  const GhostType & ghost_type,
					  const Vector<UInt> * filter_elements,
					  bool accumulate,
					  UInt id_shape,
					  UInt num_degre_of_freedom_to_interpolate,
					  UInt num_degre_of_freedom_interpolated) const {
  AKANTU_DEBUG_IN();

  Vector<Real> * shapesd_loc;
  UInt nb_element;
  UInt * conn_val;

  shapesd_loc = shapes_derivatives(type, ghost_type)[id_shape];
  nb_element  = mesh->getNbElement(type, ghost_type);
  conn_val    = mesh->getConnectivity(type, ghost_type).values;


  AKANTU_DEBUG_ASSERT(shapesd_loc != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt size_of_shapes_derivatives = ElementClass<type>::getShapeDerivativesSize();
  UInt nb_points       = control_points(type, ghost_type).getSize();
  UInt element_dimension = ElementClass<type>::getSpatialDimension();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  if(!accumulate)
    out_nablauq.clear();

  Real * shaped_val  = shapesd_loc->values;
  Real * u_val       = in_u.values;
  Real * nablauq_val = out_nablauq.values;

  UInt offset_nablauq = nb_degree_of_freedom * element_dimension;
  UInt offset_shaped  = nb_nodes_per_element * element_dimension;

  Real * shaped  = shaped_val;
  Real * u        = static_cast<Real *>(calloc(nb_nodes_per_element,
					       sizeof(Real)));

  Real * nabla_uq = static_cast<Real *>(calloc(element_dimension,
					       sizeof(Real)));


  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shaped    = shaped_val  + filter_elem_val[el] * size_of_shapes_derivatives*nb_points;
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      u[n] = u_val[conn_val[el_offset + n] * nb_degree_of_freedom + num_degre_of_freedom_to_interpolate];
    }

    for (UInt q = 0; q < nb_points; ++q) {
      /// \nabla(U) = U^t * dphi/dx
      Math::matrixt_matrix(1, element_dimension, nb_nodes_per_element,
			   u,
			   shaped,
			   nabla_uq);

      for (UInt s = 0; s < element_dimension; ++s) {
	nablauq_val[num_degre_of_freedom_interpolated * element_dimension + s] += nabla_uq[s];
      }

      nablauq_val += offset_nablauq;
      shaped      += offset_shaped;
    }
  }

  free(u);
  free(nabla_uq);

#undef INIT_VARIABLES
  AKANTU_DEBUG_OUT();
}

// /* -------------------------------------------------------------------------- */
// template <ElementType type>
// void ShapeLinked::setControlPointsByType(Vector<Real> & points){
//   control_points[type] = &points;
// }


/* -------------------------------------------------------------------------- */
/* template instanciation */
/* -------------------------------------------------------------------------- */

#define INSTANCIATE_TEMPLATE_CLASS(type)				\
  template void ShapeLinked::precomputeShapesOnControlPoints<type>(const GhostType & ghost_type); \
  template void ShapeLinked::precomputeShapeDerivativesOnControlPoints<type>(const GhostType & ghost_type); \
  template void ShapeLinked::gradientOnControlPoints<type>(const Vector<Real> &in_u, \
   							   Vector<Real> &out_nablauq, \
							   UInt num_degre_of_freedom, \
							   const GhostType & ghost_type, \
							   const Vector<UInt> * filter_elements, \
							   bool accumulate, \
							   UInt id_shape, \
							   UInt num_degre_of_freedom_to_interpolate, \
							   UInt num_degre_of_freedom_interpolated) const; \
  template void ShapeLinked::interpolateOnControlPoints<type>(const Vector<Real> &in_u, \
							      Vector<Real> &out_uq, \
							      UInt num_degre_of_freedom, \
							      const GhostType & ghost_type, \
							      const Vector<UInt> * filter_elements, \
							      bool accumulate, \
							      UInt id_shape, \
							      UInt num_degre_of_freedom_to_interpolate,	\
							      UInt num_degre_of_freedom_interpolated) const;

AKANTU_BOOST_REGULAR_ELEMENT_LIST(INSTANCIATE_TEMPLATE_CLASS)
#undef INSTANCIATE_TEMPLATE_CLASS

__END_AKANTU__
