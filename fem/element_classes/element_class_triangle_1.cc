/**
 * @file   element_class_triangle_1.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _triangle_1
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_1>::nb_nodes_per_element;
template<> UInt ElementClass<_triangle_1>::nb_quadrature_points;
template<> UInt ElementClass<_triangle_1>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_triangle_1>::shapeFunctions(const Real * x,
								 Real * shape,
								 Real * shape_deriv,
								 Real * jacobian) {
  /**
   *      t
   *      ^
   *      |
   *      x (0,0,1)
   *      |`
   *      |  `
   *      |  q `
          |  °   `
   *      x--------x-----> s
   * (1,0,0)      (0,1,0)
   *
   * N1 = 1 - s - t
   * N2 = s
   * N3 = t
   */

  Real weight = .5;

  /// shape functions
  shape[0] = 1./3.; //N1(q_0)
  shape[1] = 1./3.; //N2(q_0)
  shape[2] = 1./3.; //N3(q_0)

  /* 
   *        / dN1/ds  dN2/ds dN3/ds \
   * dnds = |                        |
   *        \ dN1/dt  dN2/dt dN3/dt /
   */
  Real dnds[nb_nodes_per_element*spatial_dimension];
  dnds[0] = -1.; dnds[1] =  1.; dnds[2] =  0.;
  dnds[3] = -1.; dnds[4] =  0.; dnds[5] =  1.;

  /// J = dxds = dnds * x
  Real dxds[spatial_dimension*spatial_dimension];
  Math::matrix_matrix(spatial_dimension, spatial_dimension, nb_nodes_per_element,
		      dnds, x, dxds);

  Real det_dxds = dxds[0] * dxds[3] - dxds[2] * dxds[1];

  /// dxds = J^{-1}
  Real inv_dxds[spatial_dimension*spatial_dimension];
  Math::inv2(dxds,inv_dxds);

  jacobian[0] = det_dxds * weight;

  Math::matrixt_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			dnds, inv_dxds, shape_deriv);
}


/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_triangle_1>::getInradius(const Real * coord) {
  return Math::triangle_inradius(coord);
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_triangle_1>::changeDimension(const Real * coord, UInt dim, UInt n_points, Real * local_coord) {
  AKANTU_DEBUG_ERROR("TO IMPLEMENT");
}
