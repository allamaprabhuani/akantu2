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
 * @section DESCRIPTION
 *
 * @verbatim
   \eta
     ^
     |
     x (0,0,1,0)
     |`
     |  `  °             \xi
     |    `    °       -
     |      `       x (0,0,0,1)
     |      q.`   -  '
     |         -`     '
     |     -       `   '
     | -             `  '
     x------------------x-----> \zeta
 (1,0,0,0)           (0,1,0,0)
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta - \zeta \\
 * N2 &=& \xi \\
 * N3 &=& \eta \\
 * N4 &=& \zeta
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f[
 * \xi_{q0} = 1/4 \qquad  \eta_{q0} = 1/4  \qquad  \zeta_{q0} = 1/4
 * @f]
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
  Real weight = 1./6.;

  /// shape functions
  shape[0] = 1./4.; /// N1(q_0)
  shape[1] = 1./4.; /// N2(q_0)
  shape[2] = 1./4.; /// N3(q_0)
  shape[3] = 1./4.; /// N4(q_0)

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccc}
   *            \frac{\partial N1}{\partial \xi}   & \frac{\partial N2}{\partial \xi}
   *          & \frac{\partial N3}{\partial \xi}   & \frac{\partial N4}{\partial \xi} \\
   *            \frac{\partial N1}{\partial \eta}  & \frac{\partial N2}{\partial \eta}
   *          & \frac{\partial N3}{\partial \eta}  & \frac{\partial N4}{\partial \eta} \\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial \zeta}
   *          & \frac{\partial N3}{\partial \zeta} & \frac{\partial N4}{\partial \zeta}
   *          \end{array}
   *        \right)
   * @f]
   */
  Real dnds[nb_nodes_per_element*spatial_dimension];
  dnds[0] = -1.; dnds[1] = 1.; dnds[2]  = 0.; dnds[3]  = 0.;
  dnds[4] = -1.; dnds[5] = 0.; dnds[6]  = 1.; dnds[7]  = 0.;
  dnds[8] = -1.; dnds[9] = 0.; dnds[10] = 0.; dnds[11] = 1.;

  /// @f$ J = dxds = dnds * x @f$
  Real dxds[spatial_dimension*spatial_dimension];
  Math::matrix_matrix(spatial_dimension, spatial_dimension, nb_nodes_per_element,
		      dnds, x, dxds);

  Real det_dxds = Math::det3(dxds);

  /// @f$ dxds = J^{-1} @f$
  Real inv_dxds[spatial_dimension*spatial_dimension];
  Math::inv3(dxds,inv_dxds);

  jacobian[0] = det_dxds * weight;
  AKANTU_DEBUG_ASSERT(jacobian[0] > 0,
		      "Negative jacobian computed, possible problem in the element node order.");

  Math::matrixt_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			dnds, inv_dxds, shape_deriv);
}


/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_tetrahedra_1>::getInradius(const Real * coord) {
  return Math::tetrahedron_inradius(coord);
}
