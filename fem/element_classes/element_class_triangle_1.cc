/**
 * @file   element_class_triangle_1.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _triangle_1
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 * @section DESCRIPTION
 *
 * @verbatim
       \eta
         ^
         |
         x (0,0,1)
         |`
         |  `
         |  q `
         |  °   `
         x--------x----->  \xi
    (1,0,0)      (0,1,0)
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta \\
 * N2 &=& \xi \\
 * N3 &=& \eta
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 1/3 \qquad  \eta_{q0} = 1/3
 * @f}
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
  Real weight = .5;

  /// shape functions
  shape[0] = 1./3.; /// N1(q_0)
  shape[1] = 1./3.; /// N2(q_0)
  shape[2] = 1./3.; /// N3(q_0)

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial \xi}  & \frac{\partial N3}{\partial \xi} \ \
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial \eta} & \frac{\partial N3}{\partial \eta}
   *          \end{array}
   *        \right)
   * @f]
   */
  Real dnds[nb_nodes_per_element*spatial_dimension];
  dnds[0] = -1.; dnds[1] =  1.; dnds[2] =  0.;
  dnds[3] = -1.; dnds[4] =  0.; dnds[5] =  1.;

  /// @f$ J = dxds = dnds * x @f$
  Real dxds[spatial_dimension*spatial_dimension];
  Math::matrix_matrix(spatial_dimension, spatial_dimension, nb_nodes_per_element,
		      dnds, x, dxds);

  Real det_dxds = dxds[0] * dxds[3] - dxds[2] * dxds[1];

  /// @f$ dxds = J^{-1} @f$
  Real inv_dxds[spatial_dimension*spatial_dimension];
  Math::inv2(dxds,inv_dxds);

  jacobian[0] = det_dxds * weight;

  AKANTU_DEBUG_ASSERT(jacobian[0] > 0,
		      "Negative jacobian computed, possible problem in the element node order.");

  Math::matrixt_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			dnds, inv_dxds, shape_deriv);
}


/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_triangle_1>::getInradius(const Real * coord) {
  return Math::triangle_inradius(coord, coord+2, coord+4);
}
