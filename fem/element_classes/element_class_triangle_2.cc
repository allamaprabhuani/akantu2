/**
 * @file   element_class_triangle_2.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _triangle_2
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
	      x 2
	      | `
	      |   `
	      |  .  `
	      |   q2  `
	    5 x          x 4
     	      |           `
     	      |             `
     	      |  .q0     q1.  `
     	      |                 `
  	      x---------x---------x-----> \xi
	      0         3         1
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{ll}
 *   \xi_{0}  = 0   &  \eta_{0} = 0   \\
 *   \xi_{1}  = 1   &  \eta_{1} = 0   \\
 *   \xi_{2}  = 0   &  \eta_{2} = 1   \\
 *   \xi_{3}  = 1/2 &  \eta_{3} = 0   \\
 *   \xi_{4}  = 1/2 &  \eta_{4} = 1/2 \\
 *   \xi_{5}  = 0   &  \eta_{5} = 1/2
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{lll}
 * N1 = -(1 - \xi - \eta) (1 - 2 (1 - \xi - \eta))
 *       & \frac{\partial N1}{\partial \xi}  = 1 - 4(1 - \xi - \eta)
 *       & \frac{\partial N1}{\partial \eta} = 1 - 4(1 - \xi - \eta) \\
 * N2 = - \xi (1 - 2 \xi)
 *       & \frac{\partial N1}{\partial \xi}  = - 1 + 4 \xi
 *       & \frac{\partial N1}{\partial \eta} = 0 \\
 * N3 = - \eta (1 - 2 \eta)
 *       & \frac{\partial N1}{\partial \xi}  = 0
 *       & \frac{\partial N1}{\partial \eta} = - 1 + 4 \eta \\
 * N4 = 4 \xi (1 - \xi - \eta)
 *       & \frac{\partial N1}{\partial \xi}  = 4 (1 - 2 \xi - \eta)
 *       & \frac{\partial N1}{\partial \eta} = - 4 \eta \\
 * N5 = 4 \xi \eta
 *       & \frac{\partial N1}{\partial \xi}  = 4 \xi
 *       & \frac{\partial N1}{\partial \eta} = 4 \eta \\
 * N6 = 4 \eta (1 - \xi - \eta)
 *       & \frac{\partial N1}{\partial \xi}  = - 4 \xi
 *       & \frac{\partial N1}{\partial \eta} = 4 (1 - \xi - 2 \eta)
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 1/6 \qquad  \eta_{q0} = 1/6 \\
 * \xi_{q1}  &=& 2/3 \qquad  \eta_{q1} = 1/6 \\
 * \xi_{q2}  &=& 1/6 \qquad  \eta_{q2} = 2/3
 * @f}
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_2>::nb_nodes_per_element;
template<> UInt ElementClass<_triangle_2>::nb_quadrature_points;
template<> UInt ElementClass<_triangle_2>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_triangle_2>::shapeFunctions(const Real * x,
								 Real * shape,
								 Real * shape_deriv,
								 Real * jacobian) {

  Real weight = 1./6.;

  /// quadrature point position
  Real quad[spatial_dimension * nb_quadrature_points];
  quad[0] = 1./6.; /// q0_{\xi}
  quad[1] = 1./6.; /// q0_{\eta}

  quad[2] = 2./3.; /// q1_{\xi}
  quad[3] = 1./6.; /// q1_{\eta}

  quad[4] = 1./6.; /// q2_{\xi}
  quad[5] = 2./3.; /// q2_{\eta}

  Real * cquad = quad;

  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    /// Natural coordinates
    Real c0 = 1 - cquad[0] - cquad[1]; /// @f$ c0 = 1 - \xi - \eta @f$
    Real c1 = cquad[0];                /// @f$ c1 = \x1 @f$
    Real c2 = cquad[1];                /// @f$ c2 = \eta @f$

    shape[0] = c0 * (2 * c0 - 1.);
    shape[1] = c1 * (2 * c1 - 1.);
    shape[2] = c2 * (2 * c2 - 1.);
    shape[3] = 4 * c0 * c1;
    shape[4] = 4 * c1 * c2;
    shape[5] = 4 * c2 * c0;

    /**
     * @f[
     * dnds =  \left(
     *   \begin{array}{cccccc}
     *       \frac{\partial N1}{\partial \xi}
     *     & \frac{\partial N2}{\partial \xi}
     *     & \frac{\partial N3}{\partial \xi}
     *     & \frac{\partial N4}{\partial \xi}
     *     & \frac{\partial N5}{\partial \xi}
     *     & \frac{\partial N6}{\partial \xi} \\
     *
     *       \frac{\partial N1}{\partial \eta}
     *     & \frac{\partial N2}{\partial \eta}
     *     & \frac{\partial N3}{\partial \eta}
     *     & \frac{\partial N4}{\partial \eta}
     *     & \frac{\partial N5}{\partial \eta}
     *     & \frac{\partial N6}{\partial \eta}
     *   \end{array}
     * \right) @f]
     */
    Real dnds[nb_nodes_per_element*spatial_dimension];
    dnds[0]  = 1 - 4 * c0;
    dnds[1]  = 4 * c1 - 1.;
    dnds[2]  = 0.;
    dnds[3]  = 4 * (c0 - c1);
    dnds[4]  = 4 * c2;
    dnds[5]  = - 4 * c2;

    dnds[6]  = 1 - 4 * c0;
    dnds[7]  = 0.;
    dnds[8]  = 4 * c2 - 1.;
    dnds[9]  = - 4 * c1;
    dnds[10] = 4 * c1;
    dnds[11] = 4 * (c0 - c2);

    /// @f$ J = dxds = dnds * x @f$
    Real dxds[spatial_dimension*spatial_dimension];
    Math::matrix_matrix(spatial_dimension, spatial_dimension, nb_nodes_per_element,
			dnds, x, dxds);

    Real det_dxds = Math::det2(dxds);

    /// @f$ dxds = J^{-1} @f$
    Real inv_dxds[spatial_dimension*spatial_dimension];
    Math::inv2(dxds, inv_dxds);

    jacobian[0] = det_dxds * weight;
    AKANTU_DEBUG_ASSERT(jacobian[0] > 0,
			"Negative jacobian computed, possible problem in the element node order.");

    Math::matrixt_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			  dnds, inv_dxds, shape_deriv);

    cquad += spatial_dimension;
    shape += nb_nodes_per_element;
    shape_deriv += nb_nodes_per_element * spatial_dimension;
    jacobian++;
  }
}


/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_triangle_2>::getInradius(const Real * coord) {
  UInt triangles[4][3] = {
    {0, 3, 5},
    {3, 1, 4},
    {3, 4, 5},
    {5, 4, 2}
  };

  Real inradius = std::numeric_limits<Real>::max();
  for (UInt t = 0; t < 4; t++) {
    Real ir = Math::triangle_inradius(coord + triangles[t][0] * 2,
				      coord + triangles[t][1] * 2,
				      coord + triangles[t][2] * 2);
    inradius = ir < inradius ? ir : inradius;
  }

  return inradius;
}
