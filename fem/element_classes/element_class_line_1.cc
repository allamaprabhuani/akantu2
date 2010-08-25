/**
 * @file   element_class_line_1.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _line_1
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

template<> UInt ElementClass<_line_1>::nb_nodes_per_element;
template<> UInt ElementClass<_line_1>::nb_quadrature_points;
template<> UInt ElementClass<_line_1>::spatial_dimension;


/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_1>::nb_nodes_per_element;
template<> UInt ElementClass<_line_1>::nb_quadrature_points;
template<> UInt ElementClass<_line_1>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_1>::shapeFunctions(const Real * x,
							     Real * shape,
							     Real * shape_deriv,
							     Real * jacobian) {
  /**
   *           q
   *  x--------|--------x
   * -1        0        1
   *
   * N1(e) = 1/2 (1 - e)
   * N2(e) = 1/2 (1 + e)
   */

  Real weight = 2;

  /// shape functions
  shape[0] = .5; //N1(0)
  shape[1] = .5; //N2(0)

  Real dnds[nb_nodes_per_element*spatial_dimension];
  ///dN1/de
  dnds[0] = - .5;
  ///dN2/de
  dnds[1] =   .5;

  Real dxds = dnds[0]*x[0] + dnds[1]*x[1];

  jacobian[0] = dxds * weight;
  AKANTU_DEBUG_ASSERT(jacobian[0]>0,
		      "negative jacobian computed problem in the element numerotation ? ");

  shape_deriv[0] = dnds[0] / dxds;
  shape_deriv[1] = dnds[1] / dxds;
}


/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_line_1>::getInradius(const Real * coord) {
  return sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_1>::changeDimension(const Real * coord, UInt dim, Real * local_coord) {
  if (dim != 2) AKANTU_DEBUG_ERROR("Is it a sens in changing coordinates from " << dim << " to line elements ?");
  Real vecs[dim*nb_nodes_per_element];
  //compute first vector
  for (UInt j = 0; j < nb_nodes_per_element-1; ++j) {
    for (UInt i = 0; i < dim; ++i) {
      vecs[i+j*dim] = coord[i+(j+1)*dim]- coord[i];
    }
  }
  //compute normal to this vector
  Math::normal2(vecs,vecs+dim);
  //normalize the vector
  Math::normalize2(vecs);
  Math::normalize2(vecs+dim);
  //invert to find the matrix rotation
  Real rotation[dim*nb_nodes_per_element];
  Math::inv2(vecs,rotation);
  //rotation of the coordinates
  //compute points
  Real points[dim*nb_nodes_per_element];
  for (UInt j = 0; j < nb_nodes_per_element; ++j) {
    for (UInt i = 0; i < dim; ++i) {
      points[i+j*dim] = coord[i+j*dim]- coord[i];
    }
  }
  Real local_coord_dim[dim*nb_nodes_per_element];
  Math::matrix_matrix(dim,dim,dim,points,rotation,local_coord_dim);
  for (UInt i = 0; i < nb_nodes_per_element; ++i) {
    local_coord[i] = local_coord_dim[i*dim];
  }
}
