/**
 * @file   shape_linked_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date   Fri Jul 15 19:41:58 2011
 *
 * @brief  ShapeLinked inline implementation
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

template <ElementKind kind>
inline void
ShapeLinked<kind>::initShapeFunctions(const Vector<Real> & nodes,
				      const types::Matrix<Real> & control_points,
				      const ElementType & type,
				      const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

#undef INIT_SHAPE_FUNCTIONS

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)					\
  setControlPointsByType<type>(control_points, ghost_type);		\
  precomputeShapesOnControlPoints<type>(nodes, ghost_type);		\
  precomputeShapeDerivativesOnControlPoints<type>(nodes, ghost_type);

template <>
inline void
ShapeLinked<_ek_structural>::initShapeFunctions(const Vector<Real> & nodes,
						const types::Matrix<Real> & control_points,
						const ElementType & type,
						const GhostType & ghost_type) {
  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

#undef INIT_SHAPE_FUNCTIONS


/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Vector<Real> & ShapeLinked<kind>::getShapes(const ElementType & type,
							 const GhostType & ghost_type,
							 UInt id) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(shapes.exists(type, ghost_type),
		      "No shapes of type "
		      << type << " in " << this->id);
  AKANTU_DEBUG_OUT();
  return *(shapes(type, ghost_type)[id]);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Vector<Real> & ShapeLinked<kind>::getShapesDerivatives(const ElementType & type,
								    const GhostType & ghost_type,
								    UInt id) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(shapes_derivatives.exists(type, ghost_type),
		      "No shapes_derivatives of type "
		      << type << " in " << this->id);
  AKANTU_DEBUG_OUT();
  return *(shapes_derivatives(type, ghost_type)[id]);
}

/* -------------------------------------------------------------------------- */
template <>
template <ElementType type>
void ShapeLinked<_ek_structural>::precomputeShapesOnControlPoints(const Vector<Real> & nodes,
								  const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();
  UInt nb_nodes_per_element           = Mesh::getNbNodesPerElement(type);

  types::Matrix<Real> & natural_coords = control_points(type, ghost_type);
  UInt nb_points = control_points(type, ghost_type).cols();

  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  std::string ghost = "";
  if(ghost_type == _ghost) {
    ghost = "ghost_";
  }

  UInt nb_element = mesh->getNbElement(type, ghost_type);
  UInt nb_shape_functions = ElementClass<type, _ek_structural>::getNbShapeFunctions();

  Vector<Real> ** shapes_tmp = new Vector<Real> *[nb_shape_functions];

  Vector<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEM::extractNodalToElementField(*mesh, nodes, x_el,
				  type, ghost_type);

  for (UInt s = 0; s < nb_shape_functions; ++s) {
    std::stringstream sstr_shapes;
    sstr_shapes << id << ":" << ghost << "shapes:" << type << ":" << s;
    shapes_tmp[s] = &(alloc<Real>(sstr_shapes.str(),
				  nb_element*nb_points,
				  size_of_shapes));
    Vector<Real>::iterator< types::Matrix<Real> > x_it = x_el.begin(spatial_dimension,
								    nb_nodes_per_element);
    Vector<Real>::iterator< types::Matrix<Real> > shapes_it =
      shapes_tmp[s]->begin_reinterpret(size_of_shapes, nb_points, nb_element);

    for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it, ++x_it) {
      types::Matrix<Real> & X = *x_it;
      types::Matrix<Real> & N = *shapes_it;
      ElementClass<type>::computeShapes(natural_coords,
					N, X,
					s);
    }
  }

  shapes(type, ghost_type) = shapes_tmp;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::precomputeShapeDerivativesOnControlPoints(const Vector<Real> & nodes,
								  const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  // Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();
  UInt natural_spatial_dimension = ElementClass<type>::getNaturalSpaceDimension();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt size_of_shapesd      = ElementClass<type>::getShapeDerivativesSize();
  types::Matrix<Real> & natural_coords = control_points(type, ghost_type);
  UInt nb_points = natural_coords.cols();

  UInt nb_element = mesh->getNbElement(type, ghost_type);
  std::string ghost = "";

  if(ghost_type == _ghost) {
    ghost = "ghost_";
  }

  Vector<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEM::extractNodalToElementField(*mesh, nodes, x_el,
				  type, ghost_type);

  UInt nb_shape_functions = ElementClass<type>::getNbShapeFunctions();

  Vector<Real> ** shapes_derivatives_tmp = new Vector<Real> *[nb_shape_functions];
  for (UInt s = 0; s < nb_shape_functions; ++s) {
    std::stringstream sstr_shapesd;
    sstr_shapesd << id << ":" << ghost << "shapes_derivatives:" << type << ":" << s;
    shapes_derivatives_tmp[s] = &(alloc<Real>(sstr_shapesd.str(),
					      nb_element*nb_points,
					      size_of_shapesd));
    Real * shapesd_val   = shapes_derivatives_tmp[s]->values;

    Vector<Real>::iterator< types::Matrix<Real> > x_it = x_el.begin(spatial_dimension,
								    nb_nodes_per_element);

    for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {
      // compute shape derivatives
      types::Matrix<Real> & X = *x_it;
      types::Tensor3<Real> B(shapesd_val,
			     natural_spatial_dimension, nb_nodes_per_element, nb_points);
      ElementClass<type>::computeShapeDerivatives(natural_coords,
						  B, X, s);

      shapesd_val += size_of_shapesd*nb_points;
    }
  }

  shapes_derivatives(type, ghost_type) = shapes_derivatives_tmp;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::extractNodalToElementField(const Vector<Real> & nodal_f,
						   Vector<Real> & elemental_f,
						   UInt num_degre_of_freedom_to_extract,
						   const GhostType & ghost_type,
						   const Vector<UInt> * filter_elements) const{
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom = nodal_f.getNbComponent();
  UInt nb_element = mesh->getNbElement(type, ghost_type);
  UInt * conn_val = mesh->getConnectivity(type, ghost_type).storage();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->storage();
  }

  elemental_f.resize(nb_element);

  Real * nodal_f_val = nodal_f.storage();
  Real * f_val = elemental_f.storage();

  UInt * el_conn;
  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) el_conn = conn_val + filter_elem_val[el] * nb_nodes_per_element;
    else el_conn = conn_val + el * nb_nodes_per_element;

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = *(el_conn + n);
      *f_val = nodal_f_val[node  * nb_degree_of_freedom + num_degre_of_freedom_to_extract];
      f_val += 1;
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::interpolateOnControlPoints(const Vector<Real> &in_u,
						   Vector<Real> &out_uq,
						   UInt nb_degree_of_freedom,
						   const GhostType & ghost_type,
						   const Vector<UInt> * filter_elements,
						   bool accumulate,
						   UInt id_shape,
						   UInt num_degre_of_freedom_to_interpolate,
						   UInt num_degre_of_freedom_interpolated) const {
  AKANTU_DEBUG_IN();

  Vector<Real> * shapes_loc = shapes(type, ghost_type)[id_shape];

  AKANTU_DEBUG_ASSERT(shapes_loc != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  Vector<Real> u_el(0, nb_nodes_per_element);
  extractNodalToElementField<type>(in_u, u_el, num_degre_of_freedom_to_interpolate,
				   ghost_type, filter_elements);

  if(!accumulate) out_uq.clear();

  UInt nb_points  = control_points(type, ghost_type).cols() * u_el.getSize();
  Vector<Real> uq(nb_points, 1, 0.);

  this->template interpolateElementalFieldOnControlPoints<type>(u_el,
								uq,
								ghost_type,
								*shapes_loc,
								filter_elements);

  for (UInt q = 0; q < nb_points; ++q) {
    out_uq(q, num_degre_of_freedom_to_interpolate) += uq(q);
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::gradientOnControlPoints(const Vector<Real> &in_u,
						Vector<Real> &out_nablauq,
						UInt nb_degree_of_freedom,
						const GhostType & ghost_type,
						const Vector<UInt> * filter_elements,
						bool accumulate,
						UInt id_shape,
						UInt num_degre_of_freedom_to_interpolate,
						UInt num_degre_of_freedom_interpolated) const {
  AKANTU_DEBUG_IN();

  Vector<Real> * shapesd_loc = shapes_derivatives(type, ghost_type)[id_shape];

  AKANTU_DEBUG_ASSERT(shapesd_loc != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  Vector<Real> u_el(0, nb_nodes_per_element);
  extractNodalToElementField<type>(in_u, u_el, num_degre_of_freedom_to_interpolate, ghost_type, filter_elements);

  UInt nb_points  = control_points(type, ghost_type).cols() * u_el.getSize();
  UInt element_dimension = ElementClass<type>::getSpatialDimension();

  Vector<Real> nablauq(nb_points, element_dimension, 0.);

  if(!accumulate) out_nablauq.clear();
  this->template gradientElementalFieldOnControlPoints<type>(u_el,
							     nablauq,
							     ghost_type,
							     *shapesd_loc,
							     filter_elements);

  Vector<Real>::iterator<types::RMatrix> nabla_u_it = nablauq.begin(1, element_dimension);
  Vector<Real>::iterator<types::RMatrix> out_nabla_u_it = out_nablauq.begin(nb_degree_of_freedom, element_dimension);
  for (UInt q = 0; q < nb_points; ++q, ++nabla_u_it, ++out_nabla_u_it) {
    for (UInt s = 0; s < element_dimension; ++s) {
      (*out_nabla_u_it)(num_degre_of_freedom_to_interpolate, s) += (*nabla_u_it)(0, s);
    }
  }

  AKANTU_DEBUG_OUT();
}
