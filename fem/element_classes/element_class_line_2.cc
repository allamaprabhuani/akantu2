/**
 * @file   element_class_line_2.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sun Oct 3 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _line_2
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 * @section DESCRIPTION
 *
 * @verbatim

             -1         0         1
  	 -----x---------x---------x-----> x
	      1         3         2
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{lll}
 *   x_{1}  = -1   &  x_{2} = 1 & x_{3} = 0
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{ll}
 *   w_1(x) = \frac{x}{2}(x - 1) & w'_1(x) = x - \frac{1}{2}\\
 *   w_2(x) = \frac{x}{2}(x + 1) & w'_2(x) = x + \frac{1}{2}\\
 *   w_3(x) =  1-x^2 & w'_3(x) = -2x
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f[
 * \begin{array}{ll}
 * x_{q1}  = -1/\sqrt{3} & x_{q2} = 1/\sqrt{3}
 * \end{array}
 * @f]
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_2>::nb_nodes_per_element;
template<> UInt ElementClass<_line_2>::nb_quadrature_points;
template<> UInt ElementClass<_line_2>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template<> inline void ElementClass<_line_2>::shapeFunctions(const Real * x,
								 Real * shape,
								 Real * shape_deriv,
								 Real * jacobian) {

  Real weight = 1;

  /// quadrature point position
  Real quad[spatial_dimension * nb_quadrature_points];
  quad[0] = -1./sqrt(3);
  quad[1] = 1./sqrt(3);

  Real * cquad = quad;

  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    shape[0] = (cquad[0] - 1) * cquad[0] / 2;
    shape[1] = (cquad[0] + 1) * cquad[0] / 2;
    shape[2] = 1 - cquad[0] * cquad[0];


    Real dnds[nb_nodes_per_element*spatial_dimension];
    dnds[0]  = cquad[0] - .5;
    dnds[1]  = cquad[0] + .5;
    dnds[2]  = -2 * cquad[0];

    Real dxds =  dnds[0] * x[0] + dnds[1] * x[1] + dnds[2] * x[2];

    jacobian[0] = dxds * weight;
    AKANTU_DEBUG_ASSERT(jacobian[0] > 0,
			"Negative jacobian computed, possible problem in the element node order.");

    shape_deriv[0] = dnds[0] / dxds;
    shape_deriv[1] = dnds[1] / dxds;
    shape_deriv[2] = dnds[2] / dxds;

    cquad += spatial_dimension;
    shape += nb_nodes_per_element;
    shape_deriv += nb_nodes_per_element * spatial_dimension;
    jacobian++;
  }
}


/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_line_2>::getInradius(const Real * coord) {
  Real dist1 = sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
  Real dist2 = sqrt((coord[1] - coord[2])*(coord[1] - coord[2]));
  return dist1 < dist2 ? dist1 : dist2;
}
