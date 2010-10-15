/**
 * @file   element_class_line_1.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _line_1
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 * @section DESCRIPTION
 *
 * @verbatim
              q
   --x--------|--------x---> x
    -1        0        1
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * w_1(x) &=& 1/2(1 - x) \\
 * w_2(x) &=& 1/2(1 + x)
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * x_{q}  &=& 0
 * @f}
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
  Real weight = 2.;

  /// shape functions
  shape[0] = .5;
  shape[1] = .5;

  Real dnds[nb_nodes_per_element*spatial_dimension];
  /// dN1/de
  dnds[0] = - .5;
  /// dN2/de
  dnds[1] =   .5;

  Real dxds = dnds[0]*x[0] + dnds[1]*x[1];

  jacobian[0] = dxds * weight;
  AKANTU_DEBUG_ASSERT(jacobian[0] > 0,
		      "Negative jacobian computed, possible problem in the element node order.");

  shape_deriv[0] = dnds[0] / dxds;
  shape_deriv[1] = dnds[1] / dxds;
}

/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_line_1>::getInradius(const Real * coord) {
  return sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_1>::translateCoordinates(const Real * coord,
								   const UInt dim,
								   const UInt n_points,
								   Real * translated_coords,
								   bool flag_inverse){
  Real coeff;
  if (flag_inverse) coeff = 1.;
  else coeff = -1.;

  for (UInt j = 0; j < n_points; ++j) {
    for (UInt i = 0; i < dim; ++i) {
      translated_coords[i+j*dim] += coeff*coord[i];
    }
  }
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_1>::computeRotationMatrix(const Real * coord,
								    const UInt dim,
								    Real * rotmatrix,
								    bool inverse_flag){
  if (dim != 2) AKANTU_DEBUG_ERROR("cannot compute non invertible rotation matrix"
				   << "=>Is it a sens in changing coordinates from "
				   << dim << " to line elements ?");
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
  if (!inverse_flag)
    Math::inv2(vecs,rotmatrix);
  else
    memcpy(rotmatrix,vecs,sizeof(Real)*dim*dim);
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_1>::changeDimension(const Real * coord,
							      const UInt dim,
							      const UInt n_points,
							      Real * local_coord) {
  Real rotation[dim*dim];
  computeRotationMatrix(coord, dim, rotation);

  Real local_coord_dim[dim*n_points];

  //rotation of the coordinates
  Real points[dim*n_points];

  memcpy(points, coord, dim * n_points * sizeof(Real));
  translateCoordinates(coord, dim, n_points, points);

  Math::matrix_matrix(dim, dim, dim, points, rotation, local_coord_dim);

  for (UInt i = 0 ; i < n_points ; ++i) {
    local_coord[i] = local_coord_dim[i*dim];
  }
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_1>::unchangeDimension(const Real * coord,
								const Real * points,
								const UInt dim,
								const UInt n_points,
								Real * local_coord) {
  Real rotation[dim*nb_nodes_per_element];
  computeRotationMatrix(coord,dim,rotation,true);

  //rotation of the coordinates
  //compute points
  Real points_dim[dim*n_points];
  memset(points_dim,0,sizeof(Real)*dim*n_points);
  for (UInt j = 0; j < n_points; ++j) {
    for (UInt i = 0; i < spatial_dimension; ++i) {
      points_dim[i+j*dim] = points[j*spatial_dimension+i];
    }
  }
  Math::matrix_matrix(n_points,dim,dim,points_dim,rotation,local_coord);
  translateCoordinates(coord,dim,n_points,local_coord,1);
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_1>::computeNormalsOnQuadPoint(const Real * coord, const UInt dim, Real * normals) {
  if (dim == 1) AKANTU_DEBUG_ERROR("cannot compute normal of line element in 1D");
  Real vec[2*dim];
  vec[0] = coord[2] - coord[0];
  vec[1] = coord[3] - coord[1];

  if (dim > 2) AKANTU_DEBUG_ERROR("To implement");
  Math::normalize2(vec);
  Math::normal2(vec,vec+dim);
  Real sign = Math::det2(vec);
  if (sign > 0){
    vec[2] *= -1;
    vec[3] *= -1;
  }

  memcpy(normals,vec+2,sizeof(Real)*2);

}
