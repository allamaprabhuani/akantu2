/**
 * @file   element_class_tetrahedra_1.cc
 * @author Guillaume ANCIAUX <anciaux@epfl.ch>
 * @date   Mon Aug 16 18:09:53 2010
 *
 * @brief  Specialization of the element_class class for the type _tetrahedra_1
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedra_1>::nb_nodes_per_element;
template<> UInt ElementClass<_tetrahedra_1>::nb_quadrature_points;
template<> UInt ElementClass<_tetrahedra_1>::spatial_dimension;

/* -------------------------------------------------------------------------- */

template<> inline void ElementClass<_tetrahedra_1>::shapeFunctions(const Real * x,
								 Real * shape,
								 Real * shape_deriv,
								 Real * jacobian) {
  /**
   *      s
   *      ^
   *      |
   *      x (0,0,1,0)
   *      |`
   *      |  `  °             t
   *      |    `    °       -  
   *      |      `       x (0,0,0,1)
   *	  |        `   -  '
   *	  |   q o   -`     '
   *	  |     -       `   '
   *	  | -             `  '
   *      x------------------x-----> r
   * (1,0,0,0)              (0,1,0,0)
   *
   * N1 = 1 - r - s - t
   * N2 = r
   * N3 = s
   * N4 = t
   */

  Real weight = .5;

  /// shape functions
  shape[0] = 1./4.; //N1(q_0)
  shape[1] = 1./4.; //N2(q_0)
  shape[2] = 1./4.; //N3(q_0)
  shape[3] = 1./4.; //N4(q_0)

  /* 
   *        / dN1/dr  dN2/dr dN3/dr dN4/dr  \
   *        |                               |
   * dnds = | dN1/ds  dN2/ds dN3/ds dN4/ds  |
   *        |                               |
   *        \ dN1/dt  dN2/dt dN3/dt dN4/dt  /
   */

  Real dnds[nb_nodes_per_element*spatial_dimension];
  dnds[0] = -1.; dnds[1] = 1.; dnds[2]  = 0.; dnds[3]  = 0.; 
  dnds[4] = -1.; dnds[5] = 0.; dnds[6]  = 1.; dnds[7]  = 0.; 
  dnds[8] = -1.; dnds[9] = 0.; dnds[10] = 0.; dnds[11] = 1.;

  /// J = dxds = dnds * x
  Real dxds[spatial_dimension*spatial_dimension];
  Math::matrix_matrix(spatial_dimension, spatial_dimension, nb_nodes_per_element,
		      dnds, x, dxds);

  Real det_dxds = Math::det3(dxds);

  /// dxds = J^{-1}
  Real inv_dxds[spatial_dimension*spatial_dimension];
  Math::inv3(dxds,inv_dxds);

  jacobian[0] = det_dxds * weight;

  Math::matrixt_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			dnds, inv_dxds, shape_deriv);
}


/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_tetrahedra_1>::getInradius(const Real * coord) {
  return Math::tetrahedron_inradius(coord);
}

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_tetrahedra_1>::changeDimension(const Real * coord, UInt dim, Real * local_coord) {
  AKANTU_DEBUG_ERROR("TO IMPLEMENT");
}



